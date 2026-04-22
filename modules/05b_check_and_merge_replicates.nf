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

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(condition), val(timepoint),
          path("*.merged.bam"),
          path("*.allMap.merged.bam"),
          path("*.spikein.merged.bam"),
          emit: merged_bams,
          optional: true

    // Passthrough: individual BAMs when not merging (concordance failed or single replicate)
    tuple val(condition), val(timepoint),
          path("passthrough/*.bam"),
          path("passthrough/*.allMap.bam"),
          path("passthrough/*.spikein.bam"),
          emit: passthrough_bams,
          optional: true

    path "concordance_report.tsv", emit: concordance_report
    path "merge_replicates.log",   emit: log

  // ── Script ────────────────────────────────────────────────────────────────
  shell:
  concordanceMin    = params.replicates?.concordance_min    ?: 0.9
  concordanceMethod = params.replicates?.concordance_method ?: 'pearson'
  sampleIdList      = sample_ids instanceof List ? sample_ids.join(' ') : sample_ids.toString()
  conditionStr      = condition.toString().replaceAll(/[^a-zA-Z0-9_-]/, '_')
  timepointStr      = timepoint.toString().replaceAll(/[^a-zA-Z0-9_-]/, '_')
  mergedPrefix      = "${conditionStr}_${timepointStr}_merged"
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

  mkdir -p passthrough

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
    for f in "${FILTERED_BAMS[@]:-}"; do [[ -f "$f" ]] && cp "$f" "passthrough/" ; done
    for f in "${ALLMAP_BAMS[@]:-}";   do [[ -f "$f" ]] && cp "$f" "passthrough/" ; done
    for f in "${SPIKE_BAMS[@]:-}";    do [[ -f "$f" ]] && cp "$f" "passthrough/" ; done
    echo -e "condition\ttimepoint\tsamples\tmin_corr\tmethod\tmerged" \
      > concordance_report.tsv
    echo -e "${CONDITION}\t${TIMEPOINT}\t${SAMPLE_IDS[*]}\tNA\t${CONCORDANCE_METHOD}\tno_single" \
      >> concordance_report.tsv
    echo "MERGE | SKIP | Single replicate — no merge needed."
    exit 0
  fi

  # ── Multi-replicate: compute concordance ───────────────────────────────────
  echo "Computing genome-wide correlation (${CONCORDANCE_METHOD})..."

  # Index BAMs if not already indexed
  for bam in "${FILTERED_BAMS[@]}"; do
    [[ -f "${bam}.bai" ]] || samtools index -@ !{task.cpus} "$bam"
  done

  # Use deepTools multiBamSummary for genome-wide binned correlation
  MULTIBAM_NPZ="multiBamSummary_bins.npz"
  CORR_TABLE="correlation_matrix.tsv"

  multiBamSummary bins \
    --bamfiles "${FILTERED_BAMS[@]}" \
    --outFileName "${MULTIBAM_NPZ}" \
    --labels "${SAMPLE_IDS[@]}" \
    --binSize 10000 \
    --numberOfProcessors !{task.cpus} \
    --extendReads 2>/dev/null || true

  # Extract correlation matrix from npz using plotCorrelation (outputs to stderr/stdout)
  plotCorrelation \
    --corData "${MULTIBAM_NPZ}" \
    --corMethod "${CONCORDANCE_METHOD}" \
    --whatToPlot scatterplot \
    --outFileCorMatrix "${CORR_TABLE}" \
    --skipZeros \
    --log1p \
    2>/dev/null || true

  # Parse minimum pairwise correlation from the matrix
  MIN_CORR="NA"
  if [[ -f "${CORR_TABLE}" ]]; then
    # Extract off-diagonal values (skip header, skip diagonal 1.0s)
    MIN_CORR=$(awk 'NR>1 {
      for (i=2; i<=NF; i++) {
        val = $i + 0
        if (val < 0.9999 && val != "") {
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

  echo -e "condition\ttimepoint\tsamples\tmin_corr\tmethod\tmerged" > concordance_report.tsv

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
      >> concordance_report.tsv
    echo "MERGE | SUCCESS | ${MERGED_PREFIX} created"

  else
    echo "✗ Concordance FAILED (${MIN_CORR} < ${CONCORDANCE_MIN}) — keeping individual replicates."
    echo "  WARNING: Replicates for ${CONDITION}@${TIMEPOINT} are too discordant to merge."
    echo "  Proceeding with individual replicates."

    for f in "${FILTERED_BAMS[@]:-}"; do [[ -f "$f" ]] && cp "$f" "passthrough/" ; done
    for f in "${ALLMAP_BAMS[@]:-}";   do [[ -f "$f" ]] && cp "$f" "passthrough/" ; done
    for f in "${SPIKE_BAMS[@]:-}";    do [[ -f "$f" ]] && cp "$f" "passthrough/" ; done

    echo -e "${CONDITION}\t${TIMEPOINT}\t${SAMPLE_IDS[*]}\t${MIN_CORR}\t${CONCORDANCE_METHOD}\tno_failed" \
      >> concordance_report.tsv
    echo "MERGE | SKIP | Concordance below threshold — individual replicates retained."
  fi
  '''
}
