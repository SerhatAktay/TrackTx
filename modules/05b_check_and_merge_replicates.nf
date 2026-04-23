// ============================================================================
// check_and_merge_replicates.nf — Replicate Concordance Check & BAM Merging
// ============================================================================
//
// Purpose:
//   Groups BAMs from the same condition+timepoint, computes genome-wide
//   pairwise correlation between replicates, and merges them if concordance
//   passes a configurable threshold.
//
// Workflow:
//   1. Group BAMs by (condition, timepoint)
//   2. For 2+ replicates: compute Pearson/Spearman correlation via
//      deepTools multiBamSummary (10kb bins) + custom awk/python correlation
//   3. If min(pairwise correlations) >= concordance_min → merge with samtools merge
//      If concordance is below threshold → warn and emit individual BAMs unchanged
//   4. Single-replicate groups are passed through unmodified
//
// Key Design Decisions:
//   • deepTools multiBamSummary bins (10kb) for fast genome-wide correlation
//   • allMap BAMs are also merged in parallel (always match filtered BAM fate)
//   • Merged sample_id is: {condition}_{timepoint}_merged
//   • Replicate metadata is collapsed (replicate set to 0 for merged samples)
//   • Spike-in BAMs are merged the same way (if present)
//
// Inputs:
//   tuple(condition, timepoint, sample_ids[], bam_files[], allmap_bams[], spike_bams[])
//
// Outputs:
//   tuple(sample_id, filtered_bam, allmap_bam, spikein_bam, condition, timepoint, replicate)
//   path  concordance_report
//
// Parameters:
//   params.replicates.concordance_min    : Minimum correlation (default: 0.9)
//   params.replicates.concordance_method : pearson | spearman (default: pearson)
//
// ============================================================================

nextflow.enable.dsl = 2

process check_and_merge_replicates {

  tag        { "${condition}_${timepoint}" }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/02_alignments/_merged",
             mode: params.publish_mode,
             overwrite: true,
             pattern: "*.{bam,bai,tsv,txt,log}"

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(condition), val(timepoint),
          val(sample_ids),
          path(bam_files),
          path(allmap_bams),
          path(spike_bams)
    val   is_paired  // "true" or "false" — controls --extendReads in multiBamSummary
    path  tss_bed   // 1-bp TSS BED (from download_genome_annotations); used to build
                    // ±500 bp TSS windows for multiBamSummary BED-file mode, which is
                    // far more appropriate for sparse PRO-seq signal than genome-wide bins

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    // Use exact filenames (not broad globs) so *.merged.bam doesn't capture
    // *.allMap.merged.bam and *.spikein.merged.bam as well, which would turn
    // filt_bam into a 3-element list and silently empty the downstream channel.
    tuple val(condition), val(timepoint),
          path("${mergedPrefix}.bam"),
          path("${mergedPrefix}.allMap.merged.bam"),
          path("${mergedPrefix}.spikein.merged.bam"),
          emit: merged_bams,
          optional: true

    // Passthrough: individual BAMs when not merging (concordance failed or single replicate).
    // Files are written to typed subdirs so the three globs are non-overlapping.
    tuple val(condition), val(timepoint),
          path("passthrough/filt/*.bam"),
          path("passthrough/allmap/*.bam"),
          path("passthrough/spike/*.bam"),
          emit: passthrough_bams,
          optional: true

    // Use condition+timepoint in the filename so each parallel invocation
    // produces a uniquely-named file.  main.nf collects and concatenates
    // these into a single cohort concordance_report.tsv via collectFile.
    // (Previously all conditions wrote to the same "concordance_report.tsv"
    // in the shared _merged/ publishDir with overwrite: true, so only the
    // last condition to finish survived.)
    path "${conditionStr}_${timepointStr}_concordance.tsv", emit: concordance_report
    path "merge_replicates.log",                            emit: log

  // ── Script ────────────────────────────────────────────────────────────────
  shell:
  concordanceMin    = params.replicates?.concordance_min    ?: 0.9
  concordanceMethod = params.replicates?.concordance_method ?: 'pearson'
  sampleIdList      = sample_ids instanceof List ? sample_ids.join(' ') : sample_ids.toString()
  conditionStr      = condition.toString().replaceAll(/[^a-zA-Z0-9_-]/, '_')
  timepointStr      = timepoint.toString().replaceAll(/[^a-zA-Z0-9_-]/, '_')
  mergedPrefix      = "${conditionStr}_${timepointStr}_merged"
  isPaired          = (is_paired?.toString() == 'true')
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  exec > >(tee -a merge_replicates.log)
  exec 2> >(tee -a merge_replicates.log >&2)

  CONDITION="!{condition}"
  TIMEPOINT="!{timepoint}"
  MERGED_PREFIX="!{mergedPrefix}"
  CONCORDANCE_MIN="!{concordanceMin}"
  CONCORDANCE_METHOD="!{concordanceMethod}"
  SAMPLE_IDS=(!{sampleIdList})

  echo "════════════════════════════════════════════════════════════"
  echo "REPLICATE MERGE | ${CONDITION} @ t=${TIMEPOINT}"
  echo "════════════════════════════════════════════════════════════"
  echo "Samples: ${SAMPLE_IDS[*]}"
  echo "Concordance method: ${CONCORDANCE_METHOD}"
  echo "Concordance threshold: ${CONCORDANCE_MIN}"
  echo ""

  # Typed passthrough subdirs so output globs don't overlap
  mkdir -p passthrough/filt passthrough/allmap passthrough/spike

  # ── Collect BAM files ──────────────────────────────────────────────────────
  # Nextflow stages files flat; reconstruct lists from staged names
  FILTERED_BAMS=()
  ALLMAP_BAMS=()
  SPIKE_BAMS=()

  for sid in "${SAMPLE_IDS[@]}"; do
    # Filtered BAMs: look for {sid}.bam (not allMap, not spikein)
    for f in *.bam; do
      [[ "$f" == *"allMap"* ]] && continue
      [[ "$f" == *"spikein"* ]] && continue
      [[ "$f" == "${sid}.bam" ]] && FILTERED_BAMS+=("$f") && break
      [[ "$f" == "${sid}_"*.bam ]] && FILTERED_BAMS+=("$f") && break
    done
    for f in *allMap*.bam; do
      [[ "$f" == *"${sid}"* ]] && ALLMAP_BAMS+=("$f") && break
    done
    for f in *spikein*.bam; do
      [[ "$f" == *"${sid}"* ]] && SPIKE_BAMS+=("$f") && break
    done
  done

  N_REPLICATES=${#FILTERED_BAMS[@]}
  echo "Detected ${N_REPLICATES} replicate BAM(s)"

  # ── Single replicate: passthrough immediately ──────────────────────────────
  if [[ "$N_REPLICATES" -le 1 ]]; then
    echo "Only 1 replicate — skipping concordance check, passing through unchanged."
    for f in "${FILTERED_BAMS[@]:-}"; do [[ -f "$f" ]] && cp "$f" "passthrough/filt/" ; done
    for f in "${ALLMAP_BAMS[@]:-}";   do [[ -f "$f" ]] && cp "$f" "passthrough/allmap/" ; done
    for f in "${SPIKE_BAMS[@]:-}";    do [[ -f "$f" ]] && cp "$f" "passthrough/spike/" ; done
    echo -e "condition\ttimepoint\tsamples\tmin_corr\tmethod\tmerged" \
      > !{conditionStr}_!{timepointStr}_concordance.tsv
    echo -e "${CONDITION}\t${TIMEPOINT}\t${SAMPLE_IDS[*]}\tNA\t${CONCORDANCE_METHOD}\tno_single" \
      >> !{conditionStr}_!{timepointStr}_concordance.tsv
    echo "MERGE | SKIP | Single replicate — no merge needed."
    exit 0
  fi

  # ── Multi-replicate: compute concordance ───────────────────────────────────
  echo "Computing genome-wide correlation (${CONCORDANCE_METHOD})..."

  # Index BAMs if not already indexed
  for bam in "${FILTERED_BAMS[@]}"; do
    [[ -f "${bam}.bai" ]] || samtools index -@ !{task.cpus} "$bam"
  done

  # ── Build ±500 bp TSS windows ─────────────────────────────────────────────
  # The 1-bp TSS BED positions are extended so that multiBamSummary BED-file
  # counts reads in the TSS ±500 bp window for each gene.
  #
  # Why BED-file mode instead of genome-wide bins?
  # PRO-seq coverage is extremely sparse (mean depth ~0.4×).  At 10 kb bins
  # virtually every bin is empty; Pearson correlation on the handful of
  # non-zero bins is meaningless (and was returning 0 for all conditions).
  # Counting directly over the 57k annotated TSS windows — where all PRO-seq
  # signal actually lives — gives a biologically meaningful correlation.
  TSS_BED="!{tss_bed}"
  TSS_WINDOWS="tss_windows.bed"
  awk 'BEGIN{OFS="\t"} NF>=3 {
    s = ($2 >= 500) ? $2 - 500 : 0
    e = $3 + 500
    print $1, s, e, (NF>=4 ? $4 : "."), (NF>=5 ? $5 : "0"), (NF>=6 ? $6 : ".")
  }' "${TSS_BED}" | LC_ALL=C sort -k1,1 -k2,2n > "${TSS_WINDOWS}"
  TSS_WIN_COUNT=$(wc -l < "${TSS_WINDOWS}" | tr -d ' ')
  echo "TSS windows built: ${TSS_WIN_COUNT} regions (±500 bp)"

  # ── deepTools: count reads per TSS window ─────────────────────────────────
  MULTIBAM_NPZ="multiBamSummary_tss.npz"
  CORR_TABLE="correlation_matrix.tsv"

  # --extendReads is only valid for paired-end libraries
  EXTEND_FLAG=""
  [[ "!{isPaired}" == "true" ]] && EXTEND_FLAG="--extendReads"

  multiBamSummary BED-file \
    --BED "${TSS_WINDOWS}" \
    --bamfiles "${FILTERED_BAMS[@]}" \
    --outFileName "${MULTIBAM_NPZ}" \
    --labels "${SAMPLE_IDS[@]}" \
    --numberOfProcessors !{task.cpus} \
    --minMappingQuality 10 \
    ${EXTEND_FLAG} 2>&1 | tee multiBamSummary.log || true

  if [[ ! -f "${MULTIBAM_NPZ}" ]]; then
    echo "ERROR: multiBamSummary did not produce ${MULTIBAM_NPZ} — check multiBamSummary.log"
    echo "Staged files in work dir:"
    ls -lh *.bam *.bai *.bed 2>/dev/null || echo "  (none found)"
  fi

  # Extract Pearson/Spearman correlation matrix.
  # Notes:
  #  • --plotFile writes a throwaway scatter PNG (required arg; /dev/null not writable here)
  #  • No --log1p            — raw counts give better Pearson on sparse PRO-seq
  #  • No --skipZeros        — keep silent TSS regions; they are informative
  plotCorrelation \
    --corData "${MULTIBAM_NPZ}" \
    --corMethod "${CONCORDANCE_METHOD}" \
    --whatToPlot scatterplot \
    --plotFile correlation_scatter.png \
    --outFileCorMatrix "${CORR_TABLE}" \
    2>&1 | tee plotCorrelation.log || true

  # Parse minimum pairwise correlation from the matrix
  MIN_CORR="NA"
  if [[ -f "${CORR_TABLE}" ]]; then
    # Extract off-diagonal values (skip header row, skip diagonal 1.0 entries)
    MIN_CORR=$(awk 'NR>1 {
      for (i=2; i<=NF; i++) {
        val = $i + 0
        if (val < 0.9999 && $i != "") {
          if (min == "" || val < min) min = val
        }
      }
    } END { print (min == "" ? "NA" : min) }' "${CORR_TABLE}")
  fi

  echo "Minimum pairwise ${CONCORDANCE_METHOD} correlation: ${MIN_CORR}"

  # ── Decide: merge or passthrough ──────────────────────────────────────────
  WILL_MERGE=false
  if [[ "${MIN_CORR}" != "NA" ]]; then
    # Use awk for floating-point comparison
    PASSES=$(awk -v corr="${MIN_CORR}" -v thr="${CONCORDANCE_MIN}" \
      'BEGIN { print (corr >= thr) ? "yes" : "no" }')
    [[ "$PASSES" == "yes" ]] && WILL_MERGE=true
  fi

  echo -e "condition\ttimepoint\tsamples\tmin_corr\tmethod\tmerged" > !{conditionStr}_!{timepointStr}_concordance.tsv

  if [[ "$WILL_MERGE" == true ]]; then
    echo "✓ Concordance passed (${MIN_CORR} >= ${CONCORDANCE_MIN}) — merging replicates."

    # Merge filtered BAMs
    samtools merge \
      -f -@ !{task.cpus} \
      "${MERGED_PREFIX}.bam" \
      "${FILTERED_BAMS[@]}"
    samtools index -@ !{task.cpus} "${MERGED_PREFIX}.bam"

    # Merge allMap BAMs
    if [[ "${#ALLMAP_BAMS[@]}" -gt 0 ]]; then
      samtools merge \
        -f -@ !{task.cpus} \
        "${MERGED_PREFIX}.allMap.merged.bam" \
        "${ALLMAP_BAMS[@]}"
      samtools index -@ !{task.cpus} "${MERGED_PREFIX}.allMap.merged.bam"
    else
      # Create empty placeholder
      touch "${MERGED_PREFIX}.allMap.merged.bam"
    fi

    # Merge spike-in BAMs (if they exist and are non-empty)
    REAL_SPIKE_BAMS=()
    for f in "${SPIKE_BAMS[@]:-}"; do
      [[ -f "$f" && -s "$f" ]] && REAL_SPIKE_BAMS+=("$f")
    done
    if [[ "${#REAL_SPIKE_BAMS[@]}" -gt 0 ]]; then
      samtools merge \
        -f -@ !{task.cpus} \
        "${MERGED_PREFIX}.spikein.merged.bam" \
        "${REAL_SPIKE_BAMS[@]}"
      samtools index -@ !{task.cpus} "${MERGED_PREFIX}.spikein.merged.bam"
    else
      touch "${MERGED_PREFIX}.spikein.merged.bam"
    fi

    echo -e "${CONDITION}\t${TIMEPOINT}\t${SAMPLE_IDS[*]}\t${MIN_CORR}\t${CONCORDANCE_METHOD}\tyes" \
      >> !{conditionStr}_!{timepointStr}_concordance.tsv
    echo "MERGE | SUCCESS | ${MERGED_PREFIX} created"

  else
    echo "✗ Concordance FAILED (${MIN_CORR} < ${CONCORDANCE_MIN}) — keeping individual replicates."
    echo "  WARNING: Replicates for ${CONDITION}@${TIMEPOINT} are too discordant to merge."
    echo "  Proceeding with individual replicates."

    for f in "${FILTERED_BAMS[@]:-}"; do [[ -f "$f" ]] && cp "$f" "passthrough/filt/" ; done
    for f in "${ALLMAP_BAMS[@]:-}";   do [[ -f "$f" ]] && cp "$f" "passthrough/allmap/" ; done
    for f in "${SPIKE_BAMS[@]:-}";    do [[ -f "$f" ]] && cp "$f" "passthrough/spike/" ; done

    echo -e "${CONDITION}\t${TIMEPOINT}\t${SAMPLE_IDS[*]}\t${MIN_CORR}\t${CONCORDANCE_METHOD}\tno_failed" \
      >> !{conditionStr}_!{timepointStr}_concordance.tsv
    echo "MERGE | SKIP | Concordance below threshold — individual replicates retained."
  fi
  '''
}
