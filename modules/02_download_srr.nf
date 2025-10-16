// ============================================================================
// download_srr.nf — SRA → FASTQ (robust, cache-aware)
// ----------------------------------------------------------------------------
// Overview
//   • (Re)use published FASTQs in ${params.output_dir}/fastq when present
//   • Best-effort prefetch → fasterq-dump (multi-thread), name normalize
//   • Optional gzip via pigz; portable checksums (md5sum|shasum)
//   • Stronger validation (non-empty, heuristic FASTQ header check)
//   • Per-sample README with provenance
//
// Workflow
//   0) Reuse cache in publishDir if files exist
//   1) prefetch (if installed; OK to fail)
//   2) fasterq-dump --split-files (-e ${task.cpus}) [--temp <dir>]
//   3) Normalize names → *_R1.fastq / *_R2.fastq
//   4) Optional: pigz compression → *.fastq.gz (+ checksums)
//   5) Sanity checks & README
//
// Inputs
//   tuple( sample_id, sra_id, condition, timepoint, replicate, paired_end )
//
// Outputs (publishDir):
//   ${params.output_dir}/01_trimmed_fastq/
//     ├── <SRR>_R1.fastq            (or .fastq.gz if compression enabled)
//     ├── <SRR>_R2.fastq            (paired only; or .fastq.gz)
//     ├── <SRR>_R1.fastq.md5|.sha256
//     ├── <SRR>_R2.fastq.md5|.sha256 (paired only)
//     └── README_fastq.txt
//
// Tunables (params.*)
//   • conda_sra     : conda spec (default: 'sra-tools>=3 pigz')
//   • advanced.sra_cpus / advanced.sra_mem
//   • sra_tmp       : tmp dir for fasterq-dump scratch (default: work dir)
//   • sra_max_size  : prefetch max size (default: '200G')
//   • fastq_gzip    : true/false (default: false)
// ============================================================================

nextflow.enable.dsl = 2

process download_srr {

  // ── Meta / resources ─────────────────────────────────────────────────────
  tag        { sra_id }
  label      'conda'
  cache      'lenient'
  // Resource allocation handled dynamically by base.config
  publishDir "${params.output_dir}/01_trimmed_fastq", mode: 'copy', overwrite: true

  // Conda env (override with --conda_sra)
  conda (params.conda_sra ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), val(sra_id),
          val(condition), val(timepoint), val(replicate),
          val(paired_end)

  // ── Declared outputs ────────────────────────────────────────────────────
  // Use glob patterns to capture either .fastq or .fastq.gz depending on params.fastq_gzip
  output:
    tuple val(sample_id),
          path("${sra_id}_R1.fastq{,.gz}"),
          path("${sra_id}_R2.fastq{,.gz}", optional: true),
          val(condition), val(timepoint), val(replicate)

    // Checksums (either md5 or sha256), optional
    path "${sra_id}*.md5",     optional: true, emit: checksums_md5
    path "${sra_id}*.sha256",  optional: true, emit: checksums_sha

    // Per-sample README (file-only)
    path "README_fastq.txt", emit: readme

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  trap 'echo "ERROR  [srr] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  export LC_ALL=C

  # ── Bindings ─────────────────────────────────────────────────────────────
  SID="!{sample_id}"
  SRR="!{sra_id}"
  COND="!{condition}"
  TP="!{timepoint}"
  REP="!{replicate}"
  PE="!{paired_end}"
  THREADS=!{task.cpus}
  OUTPUB="!{params.output_dir}/fastq"
  FQ_GZIP="!{ (params.fastq_gzip == null) ? 'false' : (params.fastq_gzip as boolean ? 'true' : 'false') }"
  SRA_TMP="!{ params.sra_tmp ?: '' }"
  SRA_MAXSZ="!{ params.sra_max_size ?: '200G' }"

  echo "INFO  [srr] ▶ sample=${SID} sra=${SRR} PE=${PE} cpus=${THREADS} gzip=${FQ_GZIP} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  mkdir -p "${OUTPUB}"

  # ── Helper: portable checksum ────────────────────────────────────────────
  mk_checksum() {
    local f="$1"
    if command -v md5sum >/dev/null 2>&1; then
      md5sum "$f" > "${f}.md5" || true
    elif command -v shasum >/dev/null 2>&1; then
      shasum -a 256 "$f" > "${f}.sha256" || true
    fi
  }

  # ── Helper: quick FASTQ validation (heuristic) ───────────────────────────
  quick_fastq_check() {
    # 1) file exists and non-empty; 2) first record header starts with '@'
    local f="$1"
    [[ -s "$f" ]] || return 1
    local head
    head=$(head -n 1 "$f" 2>/dev/null || echo "")
    [[ "$head" == @* ]] || return 2
    return 0
  }

  # ── 0) Resume from published cache (validate + link) ─────────────────────
  reuse_ok=0
  if [[ -e "${OUTPUB}/${SRR}_R1.fastq" || -e "${OUTPUB}/${SRR}_R1.fastq.gz" ]]; then
    echo "INFO  [srr] reusing published FASTQs from ${OUTPUB}"
    for cand in "${OUTPUB}/${SRR}_R1.fastq" "${OUTPUB}/${SRR}_R1.fastq.gz" \
                "${OUTPUB}/${SRR}_R2.fastq" "${OUTPUB}/${SRR}_R2.fastq.gz"; do
      if [[ -e "$cand" ]]; then
        sz=$(stat -c%s "$cand" 2>/dev/null || stat -f%z "$cand" 2>/dev/null || echo 0)
        echo "INFO  [srr]   found $(basename "$cand")  size=${sz}"
        ln -sf "$cand" .
      fi
    done
    reuse_ok=1
  fi

  # ── 1) Tools ─────────────────────────────────────────────────────────────
  command -v fasterq-dump >/dev/null 2>&1 || { echo "ERROR fasterq-dump not found" >&2; exit 1; }
  HAVE_PREFETCH=0
  if command -v prefetch >/dev/null 2>&1; then HAVE_PREFETCH=1; fi
  command -v pigz >/dev/null 2>&1 || true

  # ── 2) If no staged R1*, fetch/convert ───────────────────────────────────
  if [[ ! -e "${SRR}_R1.fastq" && ! -e "${SRR}_R1.fastq.gz" ]]; then
    # Prefetch (best-effort)
    if [[ "${HAVE_PREFETCH}" -eq 1 ]]; then
      echo "INFO  [srr] prefetch ${SRR} (max-size=${SRA_MAXSZ})"
      ( prefetch -O . --verify yes --max-size "${SRA_MAXSZ}" "${SRR}" 2>&1 | grep -v "^|" | grep -v "^2025" ) \
        || echo "WARN  [srr] prefetch failed; continuing with fasterq-dump"
    else
      echo "WARN  [srr] prefetch not available; fasterq-dump will stream"
    fi

    # Temp directory for fasterq-dump
    FQ_TMP="${SRA_TMP}"
    if [[ -n "${FQ_TMP}" ]]; then
      mkdir -p "${FQ_TMP}"
      echo "INFO  [srr] using fasterq-dump temp dir: ${FQ_TMP}"
      TMP_ARG=( --temp "${FQ_TMP}" )
    else
      TMP_ARG=()
    fi

    echo "INFO  [srr] fasterq-dump --split-files -e ${THREADS}"
    fasterq-dump --split-files -e "${THREADS}" "${TMP_ARG[@]}" -O . "${SRR}" 2>&1 | grep -E "spots read|reads read|reads written" || true

    # Normalize names (handle both *_1.fastq / *_2.fastq and single-end)
    [[ -f "${SRR}_1.fastq" ]] && mv -f "${SRR}_1.fastq" "${SRR}_R1.fastq"
    [[ -f "${SRR}_2.fastq" ]] && mv -f "${SRR}_2.fastq" "${SRR}_R2.fastq"
    [[ -f "${SRR}.fastq"    ]] && mv -f "${SRR}.fastq"    "${SRR}_R1.fastq"
  fi

  # ── 3) Optional compression with pigz ────────────────────────────────────
  if [[ "${FQ_GZIP}" == "true" ]]; then
    if command -v pigz >/dev/null 2>&1; then
      for fq in "${SRR}_R1.fastq" "${SRR}_R2.fastq"; do
        [[ -f "$fq" ]] || continue
        echo "INFO  [srr] pigz → ${fq}.gz"
        pigz -p "${THREADS}" --force "$fq"
      done
    else
      echo "WARN  [srr] pigz not found; leaving FASTQs uncompressed"
    fi
  fi

  # ── 4) Sanity checks & PE layout verification ────────────────────────────
  # Pick extensions depending on compression outcome
  R1=""
  R2=""
  if   [[ -f "${SRR}_R1.fastq.gz" ]]; then R1="${SRR}_R1.fastq.gz"
  elif [[ -f "${SRR}_R1.fastq"    ]]; then R1="${SRR}_R1.fastq"
  fi
  if   [[ -f "${SRR}_R2.fastq.gz" ]]; then R2="${SRR}_R2.fastq.gz"
  elif [[ -f "${SRR}_R2.fastq"    ]]; then R2="${SRR}_R2.fastq"
  fi

  if [[ -z "${R1}" || ! -s "${R1}" ]]; then
    echo "ERROR R1 FASTQ missing/empty for ${SRR}" >&2
    exit 1
  fi

  # Quick header check (only for uncompressed, cheap heuristic)
  if [[ "${R1}" == *.fastq ]]; then
    if ! quick_fastq_check "${R1}"; then
      echo "ERROR R1 FASTQ failed quick header check (${R1})" >&2
      exit 1
    fi
  fi

  if [[ "${PE}" == "true" ]]; then
    if [[ -z "${R2}" || ! -s "${R2}" ]]; then
      echo "ERROR Paired-end expected but R2 missing/empty for ${SRR}" >&2
      exit 1
    fi
    if [[ "${R2}" == *.fastq ]]; then
      if ! quick_fastq_check "${R2}"; then
        echo "ERROR R2 FASTQ failed quick header check (${R2})" >&2
        exit 1
      fi
    fi
  else
    # If SE run yet stray zero-byte R2 exists, drop it
    if [[ -f "${SRR}_R2.fastq" && ! -s "${SRR}_R2.fastq" ]]; then rm -f "${SRR}_R2.fastq"; fi
    if [[ -f "${SRR}_R2.fastq.gz" && ! -s "${SRR}_R2.fastq.gz" ]]; then rm -f "${SRR}_R2.fastq.gz"; fi
  fi

  # ── 5) Checksums ─────────────────────────────────────────────────────────
  mk_checksum "${R1}"
  [[ -n "${R2}" ]] && mk_checksum "${R2}"

  # ── 6) README (idempotent) ───────────────────────────────────────────────
  {
    echo "TrackTx FASTQ cache"
    echo "  sample    : ${SID}"
    echo "  accession : ${SRR}"
    echo "  layout    : ${PE}"
    echo "  gzip      : ${FQ_GZIP}"
    echo "  generated : $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
    echo
    echo "Notes"
    echo "  • Files are copied by Nextflow's publishDir for re-use on -resume."
    echo "  • If access-controlled (dbGaP/private), ensure SRA Toolkit is configured:"
    echo "      vdb-config --interactive    or    vdb-config --import <project.krt>"
  } > README_fastq.txt

  echo "INFO  [srr] ✔ done: ${R1} ${R2:+and ${R2}} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
