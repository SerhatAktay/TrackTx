// ============================================================================
// run_alignment.nf — Bowtie2 alignment (+ optional spike-in pass)
// ----------------------------------------------------------------------------
// Overview
//   • PRO-seq convention: reverse-complement R1; R2 unchanged if PE
//   • Primary genome alignment → sample_allMap.bam (sorted; all mappings)
//   • Derive primary-only BAM → sample.bam (+ .bai), QC (flagstat/idxstats)
//   • Spike-in: align residual unaligned reads (SE mode) or write empty BAM
//   • Summaries + README; unified step log
//
// Inputs
//   tuple( sample_id, read1, read2, condition, timepoint, replicate )
//   tuple( genome_id, source, genome_fa ) + path genome_bt2 (bt2 shards list)
//   tuple( spike_id,  spike_source, spike_fa ) + path spike_bt2  (bt2 shards list)
//   NOTE: SE/PE decided by params.paired_end AND presence of read2
//
// Outputs (publishDir):
//   ${params.output_dir}/02_alignments/${sample_id}/
//     • sample_allMap.bam, sample.bam/.bai, spikein.bam/.bai
//     • *.flagstat, *.idxstats, bowtie2_*.log, aligner_summary.tsv, insert_size.tsv
//     • README_alignment.txt, align_reads.log
// ============================================================================

nextflow.enable.dsl = 2

process run_alignment {

  // ── Meta / resources ─────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/02_alignments/${sample_id}", mode: 'copy', overwrite: true

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), path(read1), path(read2),
          val(condition), val(timepoint), val(replicate)

    tuple val(genome_id), val(source), path(genome_fa)
    path  genome_bt2    // MUST be the *.bt2* shards from fetch_and_build_index:index_files

    tuple val(spike_id),  val(spike_source),  path(spike_fa)
    path  spike_bt2      // optional shards; may be empty if no spike-in

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("sample.bam"),
          path("sample_allMap.bam"),
          path("spikein.bam"),
          val(condition), val(timepoint), val(replicate)

    path "bowtie2_*.log", emit: align_logs
    path "*.flagstat",    emit: flagstats
    path "*.idxstats",    emit: idxstats
    path "aligner_summary.tsv",             emit: align_summary
    path "insert_size.tsv", optional: true, emit: insert_size

    path "README_alignment.txt"
    path "align_reads.log", emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  # Unified step log
  exec > >(tee -a align_reads.log) 2>&1
  trap 'echo "ERROR  [align] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR

  ###########################################################################
  # 0) SETTINGS & SAFETY CHECKS
  ###########################################################################
  
  # Check available disk space before alignment (prevent out-of-space failures)
  AVAIL_GB=$(df -BG . 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G' || echo "0")
  if [[ "$AVAIL_GB" -lt 10 ]]; then
    echo "WARNING [align] Low disk space: ${AVAIL_GB}GB available" >&2
    echo "WARNING [align] Alignment may fail with large datasets" >&2
  fi
  THREADS=!{task.cpus}
  # More aggressive threading: use all available threads for Bowtie2
  # Leave 1 thread for SAM processing to avoid I/O bottlenecks
  BT2_T=$THREADS
  SAM_T=$(( THREADS>1 ? 1 : 1 ))

  R1="!{read1}"
  R2="!{read2}"

  # Decide SE/PE at template time from params (literal "true"/"false")
  PAIRED_PARAM='!{ params.get("paired_end", false) ? "true" : "false" }'

  IS_PE="false"
  if [[ "${PAIRED_PARAM}" == "true" && -s "${R2}" ]]; then
    IS_PE="true"
  fi

  echo "INFO  [align] ▶ sample=!{sample_id}  paired_end=${IS_PE}  cpus=${THREADS} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  if [[ "${PAIRED_PARAM}" == "true" && "${IS_PE}" != "true" ]]; then
    echo "ERROR  params.paired_end=true but R2 missing or empty" >&2; exit 1;
  fi

  ###########################################################################
  # 1) STAGE BOWTIE2 INDICES (robust to list OR single path)
  ###########################################################################
  mkdir -p bt2_index

  __BT2_GENOME_LIST__='!{ (genome_bt2 instanceof List) ? genome_bt2.collect{ it.toString() }.join(' ') : genome_bt2.toString() }'
  if [[ -n "${__BT2_GENOME_LIST__}" ]]; then
    # shellcheck disable=SC2086
    cp -f ${__BT2_GENOME_LIST__} bt2_index/
  fi
  GENOME_IDX="bt2_index/!{genome_id}"

  __BT2_SPIKE_LIST__='!{ (spike_bt2 instanceof List) ? spike_bt2.collect{ it.toString() }.join(' ') : spike_bt2.toString() }'
  SPIKE_IDX=""
  if [[ -n "!{spike_id}" && "!{spike_id}" != "none" && -n "${__BT2_SPIKE_LIST__}" ]]; then
    # shellcheck disable=SC2086
    cp -f ${__BT2_SPIKE_LIST__} bt2_index/ || true
    SPIKE_IDX="bt2_index/!{spike_id}"
  fi

  # Preflight: ensure at least one full family exists for genome
  have_full_family() {
    local p="$1"
    for ext in bt2 bt2l; do
      for s in 1 2 3 4 rev.1 rev.2; do [[ -s "${p}.${s}.${ext}" ]] || continue 2; done
      return 0
    done
    return 1
  }
  if ! have_full_family "${GENOME_IDX}"; then
    echo "ERROR  No complete Bowtie2 index shards found for prefix: ${GENOME_IDX}" >&2
    echo "DEBUG  Present files under bt2_index/:" >&2
    ls -l bt2_index/ >&2 || true
    exit 1
  fi

  ###########################################################################
  # 2) HELPERS
  ###########################################################################
  decompress() { [[ "$1" == *.gz ]] && gzip -cd -- "$1" || cat -- "$1"; }

  rc_stream() {  # reverse-complement FASTQ stream (seqkit or Python fallback)
    if command -v seqkit >/dev/null 2>&1; then
      seqkit seq --quiet -t dna -r -p -j "${THREADS}"
    else
      python3 - "$@" <<'PY'
import sys
comp=str.maketrans('ACGTNacgtn','TGCANtgcan')
it=iter(sys.stdin)
while True:
    try:
        h=next(it).rstrip('\n'); s=next(it).rstrip('\n')
        p=next(it).rstrip('\n'); q=next(it).rstrip('\n')
    except StopIteration:
        break
    sys.stdout.write(f"{h}\n{s.translate(comp)[::-1]}\n{p}\n{q[::-1]}\n")
PY
    fi
  }

  ###########################################################################
  # 3) PRIMARY GENOME ALIGNMENT  → sample_allMap.bam (sorted)
  ###########################################################################
  echo "INFO  [align] Bowtie2 → genome (!{genome_id})"

  if [[ "${IS_PE}" == "true" ]]; then
    bowtie2 -p "${BT2_T}" --end-to-end --ff --no-unal \
            -x "${GENOME_IDX}" \
            -1 <(decompress "${R2}") \
            -2 <(decompress "${R1}" | rc_stream) \
            --un-conc unaligned_R%.fastq \
            2> bowtie2_primary.log \
    | samtools sort -@ "${SAM_T}" -o sample_allMap.bam
    cat unaligned_R1.fastq unaligned_R2.fastq > unaligned.fastq
    echo "INFO  [align] primary allMap.bam ready  size=$(stat -f %z sample_allMap.bam 2>/dev/null || stat -c %s sample_allMap.bam)"
  else
    bowtie2 -p "${BT2_T}" --end-to-end --no-unal \
            -x "${GENOME_IDX}" \
            -U <(decompress "${R1}" | rc_stream) \
            --un unaligned.fastq \
            2> bowtie2_primary.log \
    | samtools sort -@ "${SAM_T}" -o sample_allMap.bam
    echo "INFO  [align] primary allMap.bam ready  size=$(stat -f %z sample_allMap.bam 2>/dev/null || stat -c %s sample_allMap.bam)"
  fi

  samtools quickcheck -v sample_allMap.bam

  ###########################################################################
  # 4) PRIMARY-ONLY BAM (keep mapped, primary; duplicates retained)
  ###########################################################################
  samtools view -@ "${THREADS}" -h -b -F 260 sample_allMap.bam \
  | samtools sort -@ "${SAM_T}" -o sample.bam

  samtools quickcheck -v sample.bam
  samtools index -@ "${THREADS}" sample.bam
  echo "INFO  [align] primary sample.bam indexed"

  samtools flagstat -@ "${THREADS}" sample.bam > sample.flagstat
  samtools idxstats        sample.bam         > sample.idxstats

  ###########################################################################
  # 5) SPIKE-IN ALIGNMENT (optional; single-end)
  ###########################################################################
  if [[ -n "${SPIKE_IDX}" && -s unaligned.fastq ]]; then
    echo "INFO  [align] Spike-in alignment (!{spike_id})"
    bowtie2 -p "${BT2_T}" --end-to-end --no-unal \
            -x "${SPIKE_IDX}" \
            -U unaligned.fastq \
            2> bowtie2_spikein.log \
    | samtools sort -@ "${SAM_T}" -o spikein.bam
    echo "INFO  [align] spikein.bam ready  size=$(stat -f %z spikein.bam 2>/dev/null || stat -c %s spikein.bam)"
  else
    echo "INFO  [align] No spike-in or no residuals → writing empty BAM"
    samtools view -@ "${THREADS}" -H sample.bam \
    | samtools view -@ "${THREADS}" -b - > spikein.bam
    : > bowtie2_spikein.log
  fi

  samtools quickcheck -v spikein.bam
  samtools index   -@ "${THREADS}" spikein.bam
  samtools flagstat -@ "${THREADS}" spikein.bam > spikein.flagstat
  samtools idxstats        spikein.bam         > spikein.idxstats

  ###########################################################################
  # 6) SUMMARIES
  ###########################################################################
  {
    echo -e "metric\tvalue"
    tot=$(awk '/in total/ {print $1}' sample.flagstat);                       [[ -n "$tot"  ]] || tot=0
    map=$(awk '/ mapped \\(/ {print $1}' sample.flagstat);                    [[ -n "$map"  ]] || map=0
    mpr=$(awk -F'[()% ]+' '/ mapped \\(/ {print $(NF-1)}' sample.flagstat);   [[ -n "$mpr"  ]] || mpr=0
    sec=$(awk '/secondary/ {print $1}' sample.flagstat);                      [[ -n "$sec"  ]] || sec=0
    dup=$(awk '/duplicates/ {print $1}' sample.flagstat);                     [[ -n "$dup"  ]] || dup=0
    echo -e "genome_total_reads\t${tot}"
    echo -e "genome_mapped_reads\t${map}"
    echo -e "genome_map_rate_pct\t${mpr}"
    echo -e "genome_secondary_reads\t${sec}"
    echo -e "genome_duplicate_reads\t${dup}"

    stot=$(awk '/in total/ {print $1}' spikein.flagstat);                     [[ -n "$stot" ]] || stot=0
    smap=$(awk '/ mapped \\(/ {print $1}' spikein.flagstat);                  [[ -n "$smap" ]] || smap=0
    smpr=$(awk -F'[()% ]+' '/ mapped \\(/ {print $(NF-1)}' spikein.flagstat); [[ -n "$smpr" ]] || smpr=0
    echo -e "spike_total_reads\t${stot}"
    echo -e "spike_mapped_reads\t${smap}"
    echo -e "spike_map_rate_pct\t${smpr}"

    if [[ "${IS_PE}" == "true" && -s bowtie2_primary.log ]]; then
      conc0=$(awk '/aligned concordantly 0 times/ {print $1}'           bowtie2_primary.log); [[ -n "$conc0" ]] || conc0=0
      conc1=$(awk '/aligned concordantly exactly 1 time/ {print $1}'    bowtie2_primary.log); [[ -n "$conc1" ]] || conc1=0
      concM=$(awk '/aligned concordantly >1 times/ {print $1}'          bowtie2_primary.log); [[ -n "$concM" ]] || concM=0
      dis=$(awk  '/aligned discordantly 1 time/ {print $1}'             bowtie2_primary.log); [[ -n "$dis"  ]] || dis=0
      echo -e "pe_concordant_0\t${conc0}"
      echo -e "pe_concordant_1\t${conc1}"
      echo -e "pe_concordant_gt1\t${concM}"
      echo -e "pe_discordant_1\t${dis}"
    fi
  } > aligner_summary.tsv

  if [[ "${IS_PE}" == "true" ]]; then
    samtools stats sample.bam \
    | awk -F '\\t' '$1=="IS"{print $2"\\t"$3}' \
    > insert_size.tsv || true
  fi

  ###########################################################################
  # 7) README
  ###########################################################################
  cat > README_alignment.txt <<TXT
Alignment artifacts for sample: !{sample_id}

Files
  • sample_allMap.bam   — all mapped alignments (primary+secondary; sorted)
  • sample.bam          — primary mapped only (-F 260); duplicates retained
  • spikein.bam         — spike-in alignments (SE); header-valid empty if none
  • sample.flagstat / sample.idxstats
  • spikein.flagstat / spikein.idxstats
  • bowtie2_primary.log / bowtie2_spikein.log
  • aligner_summary.tsv
  • insert_size.tsv (PE only)
  • align_reads.log (this step: stdout+stderr)
Notes
  - PRO-seq convention: R1 was reverse-complemented prior to alignment.
  - PE uses --ff with mate swap (-1 original R2, -2 RC(R1)).
  - Spike-in pass maps unaligned residuals in single-end mode.
TXT

  echo "INFO  [align] ✔ done ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
