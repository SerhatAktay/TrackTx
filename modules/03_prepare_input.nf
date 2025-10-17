// ============================================================================
// prepare_input.nf — Pre-processing of raw FASTQ for TrackTx-NF (single-pass)
// ----------------------------------------------------------------------------
// Overview
//   • One-pass cutadapt: adapters (-a/-A), barcodes (-u/-U), minlen pre-UMI
//   • UMI extraction (umi_tools) with 5′/3′ support
//   • Optional FastQC on raw (default: true); FastQC on final always on
//   • Emits final_R*.fastq (stub R2 for SE), QC, stats, reports, README
//
// Why it’s faster
//   • Avoids multiple cutadapt passes (adapter → barcode → pre-UMI → final)
//     by merging into a single cutadapt call before UMI extraction
//   • Skips gunzip/copy staging — tools read .gz natively
//
// Key params (with defaults):
//   params.fastqc_raw               : true   (run FastQC on input reads)
//   params.adapter_trimming.enabled : true/false
//   params.adapter_trimming.adapter1/adapter2
//   params.adapter_trimming.minlen  : user floor (optional; auto ≥ this)
//   params.barcode.enabled          : true/false
//   params.barcode.length, .location ('5'|'3'), optional R2 counterparts
//   params.umi.enabled              : true/false
//   params.umi.length, .location    : '5'|'3'
//   params.advanced.prep_insert_minlen : 12 (final insert min after UMI)
//   params.advanced.prep_cpus, prep_mem : resources
//
// Outputs (publishDir):
//   ${params.output_dir}/01_trimmed_fastq/${sample_id}/
//     ├── final_R1.fastq
//     ├── final_R2.fastq                (stub for SE to keep tuple shape)
//     ├── fastqc_raw/*                  (if enabled)
//     ├── fastqc_final/*
//     ├── trim_stats.tsv
//     ├── umi_stats.tsv                 (if UMI enabled)
//     ├── cutadapt_report.txt
//     ├── umi_extract.log               (if UMI enabled)
//     └── README_01_trimmed_fastq.txt
// ============================================================================

nextflow.enable.dsl = 2

process prepare_input {

  // ── Resources ────────────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  cache  'lenient'
  label  'conda'
  tag    { sample_id }

  // ── Clear, ordered output folder name for first-time users ──────────────
  publishDir "${params.output_dir}/01_trimmed_fastq/${sample_id}",
             mode: 'copy', overwrite: true

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(reads),                       // [R1] or [R1,R2] (gz or plain)
          val(condition), val(timepoint), val(replicate)
    val   data_type                          // "SE" or "PE"

  // ── Outputs ─────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path('final_R1.fastq'),
          path('final_R2.fastq'),            // we always create a stub for SE
          val(condition), val(timepoint), val(replicate),
          emit: cleaned_reads

    path 'fastqc_raw/*',   emit: fastqc_raw
    path 'fastqc_final/*', emit: fastqc_final

    path 'trim_stats.tsv',                emit: trim_stats
    path 'umi_stats.tsv', optional: true, emit: umi_stats

    path 'cutadapt_report.txt', emit: cutadapt_report
    path 'umi_extract.log',  optional: true, emit: umi_log

    path 'preprocess_reads.log', emit: log
    path 'README_01_trimmed_fastq.txt'

  // ── Task script (all input-derived logic lives here) ────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  trap 'echo "ERROR [prepare_input] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  export LC_ALL=C

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a preprocess_reads.log) 2>&1

  # ── Bindings from Nextflow ───────────────────────────────────────────────
  THREADS=!{ task.cpus }
  SAMPLE_ID='!{ sample_id }'
  MODE='!{ (data_type ?: "SE").toString() }'
  R1='!{ reads[0] }'
  R2='!{ (reads.size()>1 ? reads[1] : "") }'

  # Params → scalars that are safe in bash
  TRIM_ON=$([ '!{ params.adapter_trimming?.enabled == true }' = 'true' ] && echo 1 || echo 0)
  A1='!{ (params.adapter_trimming?.adapter1 ?: "").toString() }'
  A2='!{ (params.adapter_trimming?.adapter2 ?: "").toString() }'
  PRE_MINLEN_USER=!{ (params.adapter_trimming?.minlen ?: 0) }

  BC1_ON=$([ '!{ params.barcode?.enabled == true }' = 'true' ] && echo 1 || echo 0)
  BC1_LEN=!{ (params.barcode?.length  ?: 0) }
  BC1_LOC='!{ (params.barcode?.location ?: "5").toString() }'

  BC2_ON=$([ '!{ params.barcode?.r2_enabled == true }' = 'true' ] && echo 1 || echo 0)
  BC2_LEN=!{ (params.barcode?.r2_length  ?: 0) }
  BC2_LOC='!{ (params.barcode?.r2_location ?: "5").toString() }'

  UMI_ON=$([ '!{ params.umi?.enabled == true }' = 'true' ] && echo 1 || echo 0)
  UMI_LEN=!{ (params.umi?.length  ?: 0) }
  UMI_LOC='!{ (params.umi?.location ?: "5").toString() }'

  FINAL_MINLEN=!{ (params.advanced?.prep_insert_minlen ?: 12) }
  DO_FQ_RAW=$([ '!{ (params.fastqc_raw == null ? true : params.fastqc_raw) }' = 'true' ] && echo 1 || echo 0)

  # ─────────────────────────────────────────────────────────────────────────
  # Derived quantities (bash arithmetic only)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${MODE}" == "PE" && -n "${R2}" ]]; then
    BC_LOSS=$(( (BC1_ON==1 ? BC1_LEN : 0) + (BC2_ON==1 ? BC2_LEN : 0) ))
  else
    BC_LOSS=$(( (BC1_ON==1 ? BC1_LEN : 0) ))
  fi
  UMI_LOSS=$(( UMI_ON==1 ? UMI_LEN : 0 ))

  PRE_MINLEN_AUTO=$(( FINAL_MINLEN + BC_LOSS + UMI_LOSS ))
  PRE_MINLEN_EFF=$(( PRE_MINLEN_USER>PRE_MINLEN_AUTO ? PRE_MINLEN_USER : PRE_MINLEN_AUTO ))
  PRE_UMI_MIN=$(( FINAL_MINLEN + UMI_LOSS ))

  echo "INFO  [prepare_input] ▶ sample=${SAMPLE_ID} mode=${MODE} cpus=${THREADS} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "INFO  Cutadapt: trim=${TRIM_ON} adapters(A1='${A1}', A2='${A2}')"
  echo "INFO  Barcode : R1 on=${BC1_ON} len=${BC1_LEN} loc=${BC1_LOC} ; R2 on=${BC2_ON} len=${BC2_LEN} loc=${BC2_LOC}"
  echo "INFO  UMI     : on=${UMI_ON} len=${UMI_LEN} loc=${UMI_LOC}"
  echo "INFO  MinLen  : pre-UMI.effective=${PRE_MINLEN_EFF} (user=${PRE_MINLEN_USER}, auto=${PRE_MINLEN_AUTO})  final_insert_minlen=${FINAL_MINLEN}"
  QC_TOOL='!{ ((params.qc?.tool ?: "fastp").toString().toLowerCase() in ["fastp","fastqc","none"]) ? (params.qc?.tool ?: "fastp").toString().toLowerCase() : "fastp" }'
  echo "INFO  QC      : tool=${QC_TOOL} raw=${DO_FQ_RAW} final=$([[ \"${QC_TOOL}\" == \"none\" ]] && echo 0 || echo 1)"

  # ─────────────────────────────────────────────────────────────────────────
  # Temp directory handling: avoid container /tmp exhaustion
  # Use task-local writable directory for all tool temp files
  # ─────────────────────────────────────────────────────────────────────────
  TMP_BASE=".tmp"
  mkdir -p "${TMP_BASE}"
  export TMPDIR="${PWD}/${TMP_BASE}"
  export TMP="${TMPDIR}"
  export TEMP="${TMPDIR}"
  echo "INFO  TempDir : ${TMPDIR}"

  # ─────────────────────────────────────────────────────────────────────────
  # 0) Validate inputs
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${MODE}" != "SE" && "${MODE}" != "PE" ]]; then
    echo "ERROR  Mode must be SE or PE" >&2; exit 1
  fi
  if [[ "${MODE}" == "PE" && -z "${R2}" ]]; then
    echo "ERROR  Paired-end requested but R2 missing" >&2; exit 1
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # Helpers
  # ─────────────────────────────────────────────────────────────────────────
  reads_count() {  # fastq(.gz) → integer count; robust, stream-safe
    local f="$1"
    if [[ "$f" == *.gz ]]; then gzip -cd -- "$f" | awk 'END{print NR/4}'
    else                          awk 'END{print NR/4}' "$f"; fi
  }

  # ─────────────────────────────────────────────────────────────────────────
  # 1) Optional QC on RAW (fastp or FastQC)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${DO_FQ_RAW}" -eq 1 && "${QC_TOOL}" != "none" ]]; then
    mkdir -p fastqc_raw
    if [[ "${QC_TOOL}" == "fastqc" ]]; then
      echo "WARNING  Using FastQC for raw QC: slower than fastp. Consider --qc.tool fastp for speed."
      echo "INFO  FastQC (raw)"
      if [[ "${MODE}" == "PE" ]]; then
        fastqc --quiet --threads "${THREADS}" -o fastqc_raw "${R1}" "${R2}"
      else
        fastqc --quiet --threads "${THREADS}" -o fastqc_raw "${R1}"
      fi
    else
      # fastp-based QC (report only; outputs discarded)
      echo "INFO  fastp (raw QC report)"
      RAW_JSON="fastqc_raw/fastp_raw.json"
      RAW_HTML="fastqc_raw/fastp_raw.html"
      if command -v fastp >/dev/null 2>&1; then
        if [[ "${MODE}" == "PE" ]]; then
          TMP1="fastqc_raw/.raw_tmp_R1.fastq"; TMP2="fastqc_raw/.raw_tmp_R2.fastq"
          fastp -w "${THREADS}" -i "${R1}" -I "${R2}" -o "${TMP1}" -O "${TMP2}" \
                --detect_adapter_for_pe \
                --disable_adapter_trimming \
                --json "${RAW_JSON}" --html "${RAW_HTML}" --report_title "Fastp raw QC"
          rm -f "${TMP1}" "${TMP2}" || true
        else
          TMP1="fastqc_raw/.raw_tmp_R1.fastq"
          fastp -w "${THREADS}" -i "${R1}" -o "${TMP1}" \
                --disable_adapter_trimming \
                --json "${RAW_JSON}" --html "${RAW_HTML}" --report_title "Fastp raw QC"
          rm -f "${TMP1}" || true
        fi
      else
        echo "WARNING  fastp not found; falling back to FastQC for raw QC"
        if [[ "${MODE}" == "PE" ]]; then
          fastqc --quiet --threads "${THREADS}" -o fastqc_raw "${R1}" "${R2}"
        else
          fastqc --quiet --threads "${THREADS}" -o fastqc_raw "${R1}"
        fi
      fi
    fi
  elif [[ "${QC_TOOL}" == "none" ]]; then
    # Ensure directory exists for declared outputs; write a small note
    mkdir -p fastqc_raw
    echo "QC disabled (qc.tool=none). No raw QC generated." > fastqc_raw/README.txt
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 2) Single-pass cutadapt: adapters + barcodes + minlen(pre-UMI)
  #    Output → preumi_R*.fastq
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  cutadapt (one pass) → preumi_R*.fastq"
  declare -a CAD; CAD=(-j "${THREADS}")

  # Enforce pre-UMI minlen to guarantee post-UMI FINAL_MINLEN
  CAD+=(-m "${PRE_UMI_MIN}")

  # Ensure cutadapt uses the task-local temp directory
  CAD+=(--temp-dir "${TMPDIR}")

  # Adapters
  if [[ "${TRIM_ON}" -eq 1 ]]; then
    [[ -n "${A1}" ]] && CAD+=(-a "${A1}")
    if [[ "${MODE}" == "PE" && -n "${A2}" ]]; then CAD+=(-A "${A2}"); fi
  fi

  # Barcodes (signed offsets)
  if [[ "${BC1_ON}" -eq 1 && "${BC1_LEN}" -gt 0 ]]; then
    [[ "${BC1_LOC}" == "3" ]] && CAD+=(-u "-${BC1_LEN}") || CAD+=(-u "${BC1_LEN}")
  fi
  if [[ "${MODE}" == "PE" && "${BC2_ON}" -eq 1 && "${BC2_LEN}" -gt 0 ]]; then
    [[ "${BC2_LOC}" == "3" ]] && CAD+=(-U "-${BC2_LEN}") || CAD+=(-U "${BC2_LEN}")
  fi

  if [[ "${MODE}" == "PE" ]]; then
    cutadapt "${CAD[@]}" -o preumi_R1.fastq -p preumi_R2.fastq "${R1}" "${R2}" | tee cutadapt_report.txt
  else
    cutadapt "${CAD[@]}" -o preumi_R1.fastq "${R1}" | tee cutadapt_report.txt
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 3) UMI extraction → final_R*.fastq (or pass-through if UMI off)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${UMI_ON}" -eq 1 && ${UMI_LEN} -gt 0 ]]; then
    echo "INFO  UMI extraction (umi_tools extract) → final_R*.fastq"
    UMI_PAT=$(printf 'N%.0s' $(seq 1 ${UMI_LEN}))
    UMI_END_FLAG=""; [[ "${UMI_LOC}" == "3" ]] && UMI_END_FLAG="--3prime"

    if [[ "${MODE}" == "PE" ]]; then
      umi_tools extract \
        --bc-pattern="${UMI_PAT}" ${UMI_END_FLAG} \
        -I preumi_R1.fastq -S final_R1.fastq \
        --read2-in preumi_R2.fastq --read2-out final_R2.fastq \
        --log=umi_extract.log
    else
      umi_tools extract \
        --bc-pattern="${UMI_PAT}" ${UMI_END_FLAG} \
        -I preumi_R1.fastq -S final_R1.fastq \
        --log=umi_extract.log
      : > final_R2.fastq   # keep tuple shape consistent for SE
    fi
  else
    echo "INFO  UMI extraction: OFF → pass-through preumi_* to final_*"
    if [[ "${MODE}" == "PE" ]]; then
      cp -f preumi_R1.fastq final_R1.fastq
      cp -f preumi_R2.fastq final_R2.fastq
    else
      cp -f preumi_R1.fastq final_R1.fastq
      : > final_R2.fastq   # SE stub
    fi
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 4) QC on FINAL reads (fastp or FastQC)
  # ─────────────────────────────────────────────────────────────────────────
  mkdir -p fastqc_final
  if [[ "${QC_TOOL}" == "none" ]]; then
    echo "INFO  QC disabled on final reads (qc.tool=none)"
    echo "QC disabled (qc.tool=none). No final QC generated." > fastqc_final/README.txt
  elif [[ "${QC_TOOL}" == "fastqc" ]]; then
    echo "WARNING  Using FastQC for final QC: slower than fastp. Consider --qc.tool fastp for speed."
    echo "INFO  FastQC (final)"
    if [[ "${MODE}" == "PE" ]]; then
      fastqc --quiet --threads "${THREADS}" -o fastqc_final final_R1.fastq final_R2.fastq
    else
      fastqc --quiet --threads "${THREADS}" -o fastqc_final final_R1.fastq
    fi
  else
    echo "INFO  fastp (final QC report)"
    FINAL_JSON="fastqc_final/fastp_final.json"
    FINAL_HTML="fastqc_final/fastp_final.html"
    if command -v fastp >/dev/null 2>&1; then
      if [[ "${MODE}" == "PE" ]]; then
        TMP1="fastqc_final/.final_tmp_R1.fastq"; TMP2="fastqc_final/.final_tmp_R2.fastq"
        fastp -w "${THREADS}" -i final_R1.fastq -I final_R2.fastq -o "${TMP1}" -O "${TMP2}" \
              --detect_adapter_for_pe \
              --disable_adapter_trimming \
              --json "${FINAL_JSON}" --html "${FINAL_HTML}" --report_title "Fastp final QC"
        rm -f "${TMP1}" "${TMP2}" || true
      else
        TMP1="fastqc_final/.final_tmp_R1.fastq"
        fastp -w "${THREADS}" -i final_R1.fastq -o "${TMP1}" \
              --disable_adapter_trimming \
              --json "${FINAL_JSON}" --html "${FINAL_HTML}" --report_title "Fastp final QC"
        rm -f "${TMP1}" || true
      fi
    else
      echo "WARNING  fastp not found; falling back to FastQC for final QC"
      if [[ "${MODE}" == "PE" ]]; then
        fastqc --quiet --threads "${THREADS}" -o fastqc_final final_R1.fastq final_R2.fastq
      else
        fastqc --quiet --threads "${THREADS}" -o fastqc_final final_R1.fastq
      fi
    fi
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 5) Stats: trim_stats.tsv (input → preUMI), umi_stats.tsv (preUMI → final)
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  Computing read-count stats"
  R1_IN=$(reads_count "${R1}")
  R2_IN=$([[ "${MODE}" == "PE" ]] && reads_count "${R2}" || echo 0)

  R1_PRE=$(reads_count preumi_R1.fastq)
  R2_PRE=$([[ "${MODE}" == "PE" ]] && reads_count preumi_R2.fastq || echo 0)

  R1_FIN=$(reads_count final_R1.fastq)
  R2_FIN=$([[ "${MODE}" == "PE" ]] && reads_count final_R2.fastq || echo 0)

  {
    echo -e "sample_id\tread\treads_in\treads_preumi\tpct_kept_preumi"
    printf "%s\tR1\t%s\t%s\t%.4f\n" "${SAMPLE_ID}" "${R1_IN}" "${R1_PRE}" $(awk -v a="${R1_IN}" -v b="${R1_PRE}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
    if [[ "${MODE}" == "PE" ]]; then
      printf "%s\tR2\t%s\t%s\t%.4f\n" "${SAMPLE_ID}" "${R2_IN}" "${R2_PRE}" $(awk -v a="${R2_IN}" -v b="${R2_PRE}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
    fi
  } > trim_stats.tsv

  if [[ "${UMI_ON}" -eq 1 && ${UMI_LEN} -gt 0 ]]; then
    {
      echo -e "sample_id\tumi_on\tumi_len\tumi_loc\tread\treads_preumi\treads_final\tpct_kept_final"
      printf "%s\t1\t%s\t%s\tR1\t%s\t%s\t%.4f\n" "${SAMPLE_ID}" "${UMI_LEN}" "${UMI_LOC}" "${R1_PRE}" "${R1_FIN}" $(awk -v a="${R1_PRE}" -v b="${R1_FIN}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
      if [[ "${MODE}" == "PE" ]]; then
        printf "%s\t1\t%s\t%s\tR2\t%s\t%s\t%.4f\n" "${SAMPLE_ID}" "${UMI_LEN}" "${UMI_LOC}" "${R2_PRE}" "${R2_FIN}" $(awk -v a="${R2_PRE}" -v b="${R2_FIN}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
      fi
    } > umi_stats.tsv
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 6) README
  # ─────────────────────────────────────────────────────────────────────────
  cat > README_01_trimmed_fastq.txt <<TXT
Prepared (trimmed/UMI-processed) reads for sample: !{sample_id}

Files
  - final_R1.fastq            : cleaned R1 (post cutadapt + UMI)
  - final_R2.fastq            : cleaned R2 (stub for SE; real for PE)
  - fastqc_raw/               : FastQC on raw inputs (if enabled)
  - fastqc_final/             : FastQC on FINAL reads
  - trim_stats.tsv            : read counts (input → preUMI)
  - umi_stats.tsv             : read counts (preUMI → final; if UMI enabled)
  - cutadapt_report.txt       : full cutadapt report
  - umi_extract.log           : umi_tools log (if UMI enabled)
  - preprocess_reads.log      : step log (stdout+stderr)

Parameters
  fastqc_raw                = ${DO_FQ_RAW}
  adapter_trim.enabled      = ${TRIM_ON} (A1='${A1}' A2='${A2}')
  adapter_trim.minlen(user) = ${PRE_MINLEN_USER}
  preUMI.minlen(effective)  = ${PRE_UMI_MIN}
  barcode.R1                = ${BC1_ON} (len=${BC1_LEN} loc=${BC1_LOC})
  barcode.R2                = ${BC2_ON} (len=${BC2_LEN} loc=${BC2_LOC})
  umi                       = ${UMI_ON} (len=${UMI_LEN} loc=${UMI_LOC})
  final_insert_minlen       = ${FINAL_MINLEN}
TXT

  echo "INFO  [prepare_input] ✔ done ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
