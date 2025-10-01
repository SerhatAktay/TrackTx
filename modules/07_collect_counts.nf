// ============================================================================
// 07_collect_counts.nf — Total mapped reads (main, allMap, spike-in)
// ----------------------------------------------------------------------------
// Overview
//   • Counts mapped reads from:
//       - MAIN (filtered primary BAM)
//       - ALLMAP (all mapped alignments BAM; 0 if missing/empty)
//       - SPIKE (optional spike-in BAM or '-' literal → 0)
//   • Ensures BAMs are indexed (configurable) and writes a one-row TSV.
//
// Why it’s faster
//   • Skips index building by default (fail fast if .bai missing).
//   • Uses idxstats (constant-time per contig; minimal I/O) once per BAM.
//
// Inputs
//   tuple( sid, main_bam, allmap_bam, spike_bam_or_dash, condition, timepoint, replicate )
//
// Outputs (publishDir):
//   ${params.output_dir}/04_counts/${sid}/
//     ├── ${sid}.counts.tsv
//     └── collect_counts.log
//
// Config knobs
//   • params.counts_allow_index_build : true (default in nextflow.config)
//       - false → require .bai to exist; error if missing
//       - true  → build index on the fly when missing
// ============================================================================

nextflow.enable.dsl = 2

process collect_counts {

  // ── Meta / resources ─────────────────────────────────────────────────────
  tag        { sid }
  label      'conda'
  cache      'lenient'
  // Resource allocation handled dynamically by base.config

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/04_counts/${sid}", mode: 'copy', overwrite: true

  // Optional conda (samtools only)
  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sid), path(main_bam), path(allmap_bam), val(spike_in), val(cond), val(tp), val(rep)

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    path "${sid}.counts.tsv"
    path "collect_counts.log", optional: true, emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  trap 'echo "ERROR  [counts] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a collect_counts.log) 2>&1

  # ── Bindings / locals ────────────────────────────────────────────────────
  SID='!{sid}'
  MAIN_BAM='!{main_bam}'
  ALLMAP_BAM='!{allmap_bam}'
  SPIKE_IN='!{spike_in}'
  THREADS=!{task.cpus}
  ALLOW_INDEX_STR='!{ params.get("counts_allow_index_build", false) ? "true" : "false" }'

  # Convert to 1/0 without $(...)
  if [[ "${ALLOW_INDEX_STR}" == "true" ]]; then ALLOW_INDEX=1; else ALLOW_INDEX=0; fi

  echo "INFO  [counts] ▶ sample=${SID}  allow_index_build=${ALLOW_INDEX} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  command -v samtools >/dev/null 2>&1 || { echo "ERROR samtools not found" >&2; exit 1; }

  # ─────────────────────────────────────────────────────────────────────────
  # Helpers
  # ─────────────────────────────────────────────────────────────────────────
  require_or_index () {
    # $1=bam
    local bam="$1"
    [[ -s "$bam" ]] || return 1
    local bai="${bam}.bai"
    if [[ -s "$bai" ]]; then
      return 0
    fi
    if [[ "$ALLOW_INDEX" -eq 1 ]]; then
      echo "INFO  indexing: ${bam}"
      samtools index -@ "${THREADS}" "$bam"
      [[ -s "$bai" ]] || { echo "ERROR failed to create index for ${bam}" >&2; return 1; }
    else
      echo "ERROR  missing index: ${bam}.bai (set params.counts_allow_index_build=true to auto-build)" >&2
      return 1
    fi
  }

  sum_idxstats () {
    # $1=bam → prints total mapped across contigs (excludes '*')
    local bam="$1"
    samtools idxstats "$bam" | awk '$1!="*" {m+=$3} END{print (m?m:0)}'
  }

  # ─────────────────────────────────────────────────────────────────────────
  # 1) MAIN (filtered) mapped reads
  # ─────────────────────────────────────────────────────────────────────────
  [[ -s "${MAIN_BAM}" ]] || { echo "ERROR main BAM missing/empty: ${MAIN_BAM}" >&2; exit 1; }
  require_or_index "${MAIN_BAM}"
  MAIN=$(sum_idxstats "${MAIN_BAM}")

  # ─────────────────────────────────────────────────────────────────────────
  # 2) ALLMAP (all mapped) reads
  # ─────────────────────────────────────────────────────────────────────────
  if [[ -s "${ALLMAP_BAM}" ]]; then
    require_or_index "${ALLMAP_BAM}"
    ALLMAP=$(sum_idxstats "${ALLMAP_BAM}")
  else
    echo "WARN  allMap BAM missing/empty → allmap_reads=0"
    ALLMAP=0
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 3) SPIKE-IN mapped reads (optional)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${SPIKE_IN}" != "-" && -s "${SPIKE_IN}" ]]; then
    # header-only empties still yield 0 via idxstats; index if allowed
    if [[ -s "${SPIKE_IN}.bai" || "${ALLOW_INDEX}" -eq 1 ]]; then
      require_or_index "${SPIKE_IN}" || true
    fi
    SPIK=$(sum_idxstats "${SPIKE_IN}")
  else
    SPIK=0
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 4) Write TSV (header + one row)
  # ─────────────────────────────────────────────────────────────────────────
  {
    printf "sample\tmain_reads\tallmap_reads\tspike_reads\treplicate\tcondition\ttimepoint\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${SID}" "${MAIN}" "${ALLMAP}" "${SPIK}" "!{rep}" "!{cond}" "!{tp}"
  } > "${SID}.counts.tsv"

  echo "INFO  [counts] ${SID}  main=${MAIN}  allMap=${ALLMAP}  spike=${SPIK} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
