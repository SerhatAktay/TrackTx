// ============================================================================
// prepare_input.nf — FASTQ Preprocessing and Quality Control
// ============================================================================
//
// Purpose:
//   Single-pass preprocessing of raw FASTQ files for PRO-seq analysis
//
// Features:
//   • One-pass cutadapt: adapters, barcodes, and length filtering
//   • UMI extraction with 5'/3' support (umi_tools)
//   • Quality control on raw and final reads (FastQC)
//   • Comprehensive statistics and reporting
//   • Fast execution with multi-threading
//
// Why Single-Pass?
//   Traditional multi-pass approach:
//     raw → adapter trim → barcode trim → UMI extract → final
//   
//   Our optimized approach:
//     raw → cutadapt (adapters + barcodes + minlen) → UMI extract → final
//   
//   Benefits: Faster execution, fewer intermediate files, simpler logic
//
// Workflow:
//   1. Optional QC on raw reads (FastQC)
//   2. Single-pass cutadapt (adapters + barcodes + length filter)
//   3. UMI extraction (if enabled)
//   4. QC on final cleaned reads (FastQC)
//   5. Generate statistics and reports
//
// Inputs:
//   tuple(sample_id, reads, condition, timepoint, replicate)
//   val(data_type) : "SE" or "PE"
//
// Outputs:
//   ${params.output_dir}/01_trimmed_fastq/${sample_id}/
//     ├── final_R1.fastq              — Cleaned R1
//     ├── final_R2.fastq              — Cleaned R2 (or stub for SE)
//     ├── fastqc_raw/                 — Raw QC reports
//     ├── fastqc_final/               — Final QC reports
//     ├── trim_stats.tsv              — Trimming statistics
//     ├── umi_stats.tsv               — UMI extraction stats (if enabled)
//     ├── cutadapt_report.txt         — Detailed cutadapt log
//     ├── umi_extract.log             — UMI extraction log (if enabled)
//     ├── preprocess_reads.log        — Complete process log
//     └── README_01_trimmed_fastq.txt — Documentation
//
// Parameters:
//   params.fastqc_raw                     : Run QC on raw reads (default: true)
//   params.qc.enabled                     : Enable QC (default: true)
//   params.adapter_trimming.enabled       : Enable adapter trimming
//   params.adapter_trimming.adapter1/2    : Adapter sequences
//   params.adapter_trimming.minlen        : Minimum length (user override)
//   params.barcode.enabled                : Enable barcode removal
//   params.barcode.length/location        : Barcode specs (5' or 3')
//   params.umi.enabled                    : Enable UMI extraction
//   params.umi.length/location            : UMI specs (5' or 3')
//   params.advanced.prep_insert_minlen    : Final insert minimum (default: 12)
//
// ============================================================================

nextflow.enable.dsl = 2

process prepare_input {

  tag    { sample_id }
  label  'conda'
  cache  { params.prepare_input_lenient_cache ? 'lenient' : 'deep' }

  publishDir "${params.output_dir}/01_trimmed_fastq/${sample_id}",
             mode: 'copy',
             overwrite: true,
             saveAs: { filename ->
               def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
               // Only publish processed outputs, NOT the raw input FASTQ files
               // This prevents duplicating large raw FASTQ files in the results folder
               // Publish only: final_*.fastq, QC reports, stats, logs, and README
               if (name.startsWith('final_') || 
                   name.startsWith('fastqc_') || 
                   name.contains('_stats.tsv') ||
                   name.contains('_report.txt') ||
                   name.contains('_extract.log') ||
                   name.contains('preprocess_reads.log') ||
                   name.startsWith('README_')) {
                 return name
               }
               // Skip everything else (raw input FASTQ files)
               return null
             }

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(reads),
          val(condition), val(timepoint), val(replicate)
    val   data_type

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path('final_R1.fastq'),
          path('final_R2.fastq'),
          val(condition), val(timepoint), val(replicate),
          emit: cleaned_reads

    path 'fastqc_raw/*',                  emit: fastqc_raw
    path 'fastqc_final/*',                emit: fastqc_final
    path 'trim_stats.tsv',                emit: trim_stats
    path 'umi_stats.tsv',    optional: true, emit: umi_stats
    path 'cutadapt_report.txt',           emit: cutadapt_report
    path 'umi_extract.log',  optional: true, emit: umi_log
    path 'preprocess_reads.log',          emit: log
    path 'README_01_trimmed_fastq.txt'

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a preprocess_reads.log) 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "PREP | START | sample=!{sample_id} | mode=!{data_type ?: "SE"} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID='!{sample_id}'
  THREADS=!{task.cpus}
  MODE='!{(data_type ?: "SE").toString()}'
  
  R1='!{reads[0]}'
  R2='!{(reads.size() > 1 ? reads[1] : "")}'

  echo "PREP | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "PREP | CONFIG | Mode: ${MODE}"
  echo "PREP | CONFIG | Threads: ${THREADS}"
  echo "PREP | CONFIG | R1: ${R1}"
  if [[ -n "${R2}" ]]; then
    echo "PREP | CONFIG | R2: ${R2}"
  fi

  # ── Adapter Trimming Parameters ──
  TRIM_ENABLED=$([ '!{params.adapter_trimming?.enabled == true}' = 'true' ] && echo 1 || echo 0)
  ADAPTER1='!{(params.adapter_trimming?.adapter1 ?: "").toString()}'
  ADAPTER2='!{(params.adapter_trimming?.adapter2 ?: "").toString()}'
  USER_MINLEN=!{(params.adapter_trimming?.minlen ?: 0)}

  echo "PREP | CONFIG | Adapter trimming: $([ ${TRIM_ENABLED} -eq 1 ] && echo "enabled" || echo "disabled")"
  if [[ ${TRIM_ENABLED} -eq 1 ]]; then
    echo "PREP | CONFIG | Adapter 1: ${ADAPTER1:-none}"
    echo "PREP | CONFIG | Adapter 2: ${ADAPTER2:-none}"
  fi

  # ── Barcode Parameters ──
  BC1_ENABLED=$([ '!{params.barcode?.enabled == true}' = 'true' ] && echo 1 || echo 0)
  BC1_LENGTH=!{(params.barcode?.length ?: 0)}
  BC1_LOCATION='!{(params.barcode?.location ?: "5").toString()}'

  BC2_ENABLED=$([ '!{params.barcode?.r2_enabled == true}' = 'true' ] && echo 1 || echo 0)
  BC2_LENGTH=!{(params.barcode?.r2_length ?: 0)}
  BC2_LOCATION='!{(params.barcode?.r2_location ?: "5").toString()}'

  echo "PREP | CONFIG | Barcode R1: $([ ${BC1_ENABLED} -eq 1 ] && echo "enabled (${BC1_LENGTH}bp, ${BC1_LOCATION}')" || echo "disabled")"
  if [[ "${MODE}" == "PE" ]]; then
    echo "PREP | CONFIG | Barcode R2: $([ ${BC2_ENABLED} -eq 1 ] && echo "enabled (${BC2_LENGTH}bp, ${BC2_LOCATION}')" || echo "disabled")"
  fi

  # ── UMI Parameters ──
  UMI_ENABLED=$([ '!{params.umi?.enabled == true}' = 'true' ] && echo 1 || echo 0)
  UMI_LENGTH=!{(params.umi?.length ?: 0)}
  UMI_LOCATION='!{(params.umi?.location ?: "5").toString()}'

  echo "PREP | CONFIG | UMI: $([ ${UMI_ENABLED} -eq 1 ] && echo "enabled (${UMI_LENGTH}bp, ${UMI_LOCATION}')" || echo "disabled")"

  # ── QC Parameters ──
  QC_ENABLED=$([ '!{(params.qc?.enabled == null ? true : params.qc?.enabled)}' = 'true' ] && echo 1 || echo 0)
  QC_RAW=$([ '!{(params.fastqc_raw == null ? true : params.fastqc_raw)}' = 'true' ] && echo 1 || echo 0)
  FINAL_MINLEN=!{(params.advanced?.prep_insert_minlen ?: 12)}

  echo "PREP | CONFIG | QC (FastQC): $([ ${QC_ENABLED} -eq 1 ] && echo "enabled" || echo "disabled")"
  echo "PREP | CONFIG | QC on raw reads: $([ ${QC_RAW} -eq 1 ] && echo "yes" || echo "no")"
  echo "PREP | CONFIG | Final minimum length: ${FINAL_MINLEN}bp"

  # ── Calculate Length Requirements ──
  if [[ "${MODE}" == "PE" && -n "${R2}" ]]; then
    BARCODE_LOSS=$(( (BC1_ENABLED == 1 ? BC1_LENGTH : 0) + (BC2_ENABLED == 1 ? BC2_LENGTH : 0) ))
  else
    BARCODE_LOSS=$(( (BC1_ENABLED == 1 ? BC1_LENGTH : 0) ))
  fi
  
  UMI_LOSS=$(( UMI_ENABLED == 1 ? UMI_LENGTH : 0 ))
  AUTO_MINLEN=$(( FINAL_MINLEN + BARCODE_LOSS + UMI_LOSS ))
  EFFECTIVE_MINLEN=$(( USER_MINLEN > AUTO_MINLEN ? USER_MINLEN : AUTO_MINLEN ))
  PRE_UMI_MINLEN=$(( FINAL_MINLEN + UMI_LOSS ))

  echo "PREP | CONFIG | Length calculations:"
  echo "PREP | CONFIG |   Barcode loss: ${BARCODE_LOSS}bp"
  echo "PREP | CONFIG |   UMI loss: ${UMI_LOSS}bp"
  echo "PREP | CONFIG |   Auto min length (pre-UMI): ${AUTO_MINLEN}bp"
  echo "PREP | CONFIG |   User min length: ${USER_MINLEN}bp"
  echo "PREP | CONFIG |   Effective min length (pre-UMI): ${EFFECTIVE_MINLEN}bp"
  echo "PREP | CONFIG |   Min length (post-UMI): ${PRE_UMI_MINLEN}bp"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "PREP | VALIDATE | Checking inputs..."

  # Fix: Files may be gzipped but have .fastq extension.
  # Case 1: Symlink from download_srr (SRR_R1.fastq -> SRR_R1.fastq.gz). Use target, do NOT mv
  #   (mv would overwrite the real .gz file with the symlink and break it).
  # Case 2: Regular file with .fastq extension but gzip magic bytes. Rename to .gz.
  fix_gzip_extension() {
    local f="$1"
    [[ -z "$f" || ! -e "$f" ]] && echo "$f" && return
    [[ "$f" == *.gz ]] && echo "$f" && return
    # Symlink: use target (the actual .gz file) - do not mv
    if [[ -L "$f" ]]; then
      local dir base target
      dir=$(dirname "$f")
      base=$(basename "$f")
      target=$(readlink "$f")  # target may be relative
      if [[ "$target" != /* ]]; then
        target="${dir}/${target}"
      fi
      if [[ -f "$target" && "$target" == *.gz ]]; then
        echo "PREP | VALIDATE | Detected symlink .fastq -> .gz, using target: $target"
        echo "$target"
        return
      fi
    fi
    # Regular file: check gzip magic bytes
    local magic
    magic=$(head -c 2 "$f" 2>/dev/null | od -A n -t x1 2>/dev/null | tr -d ' \n' | head -c 4)
    if [[ "$magic" == "1f8b" ]]; then
      echo "PREP | VALIDATE | Detected gzipped content with .fastq extension, renaming to .gz: $f"
      mv "$f" "${f}.gz"
      echo "${f}.gz"
    else
      echo "$f"
    fi
  }
  R1_FIXED=$(fix_gzip_extension "${R1}")
  R1="${R1_FIXED}"
  if [[ -n "${R2}" && -f "${R2}" ]]; then
    R2_FIXED=$(fix_gzip_extension "${R2}")
    R2="${R2_FIXED}"
  fi

  if [[ "${MODE}" != "SE" && "${MODE}" != "PE" ]]; then
    echo "PREP | ERROR | Mode must be SE or PE, got: ${MODE}"
    exit 1
  fi

  if [[ "${MODE}" == "PE" && -z "${R2}" ]]; then
    echo "PREP | ERROR | Paired-end mode but R2 file is missing"
    exit 1
  fi

  if [[ ! -f "${R1}" ]]; then
    echo "PREP | ERROR | R1 file not found: ${R1}"
    exit 1
  fi

  if [[ "${MODE}" == "PE" && ! -f "${R2}" ]]; then
    echo "PREP | ERROR | R2 file not found: ${R2}"
    exit 1
  fi

  # Check file sizes
  R1_SIZE=$(stat -c%s "${R1}" 2>/dev/null || stat -f%z "${R1}" 2>/dev/null || echo "unknown")
  echo "PREP | VALIDATE | R1 size: ${R1_SIZE} bytes"
  
  if [[ "${MODE}" == "PE" ]]; then
    R2_SIZE=$(stat -c%s "${R2}" 2>/dev/null || stat -f%z "${R2}" 2>/dev/null || echo "unknown")
    echo "PREP | VALIDATE | R2 size: ${R2_SIZE} bytes"
  fi

  echo "PREP | VALIDATE | Input validation complete"

  ###########################################################################
  # 3) HELPER FUNCTIONS
  ###########################################################################

  # Count reads in FASTQ file (handles .gz)
  count_reads() {
    local file="$1"
    if [[ "$file" == *.gz ]]; then
      gzip -cd "$file" | awk 'END{print NR/4}' 2>/dev/null || echo 0
    else
      awk 'END{print NR/4}' "$file" 2>/dev/null || echo 0
    fi
  }

  ###########################################################################
  # 4) QC ON RAW READS (Optional)
  ###########################################################################

  mkdir -p fastqc_raw

  if [[ ${QC_RAW} -eq 1 && ${QC_ENABLED} -eq 1 ]]; then
    echo "PREP | QC-RAW | Running FastQC on raw reads..."
    
    if [[ "${MODE}" == "PE" ]]; then
      fastqc --quiet --threads "${THREADS}" -o fastqc_raw "${R1}" "${R2}"
    else
      fastqc --quiet --threads "${THREADS}" -o fastqc_raw "${R1}"
    fi
    
    echo "PREP | QC-RAW | FastQC reports generated"
  else
    echo "PREP | QC-RAW | Skipped (disabled)"
    echo "QC on raw reads disabled." > fastqc_raw/README.txt
  fi

  ###########################################################################
  # 5) SINGLE-PASS CUTADAPT
  ###########################################################################

  echo "PREP | CUTADAPT | Running single-pass trimming..."
  echo "PREP | CUTADAPT | Output: preumi_R*.fastq"

  # Build cutadapt command
  CUTADAPT_CMD=("-j" "${THREADS}")
  
  # Minimum length filter (pre-UMI)
  CUTADAPT_CMD+=("-m" "${PRE_UMI_MINLEN}")
  echo "PREP | CUTADAPT | Minimum length filter: ${PRE_UMI_MINLEN}bp"

  # Adapter trimming
  if [[ ${TRIM_ENABLED} -eq 1 ]]; then
    if [[ -n "${ADAPTER1}" ]]; then
      CUTADAPT_CMD+=("-a" "${ADAPTER1}")
      echo "PREP | CUTADAPT | R1 adapter: ${ADAPTER1}"
    fi
    if [[ "${MODE}" == "PE" && -n "${ADAPTER2}" ]]; then
      CUTADAPT_CMD+=("-A" "${ADAPTER2}")
      echo "PREP | CUTADAPT | R2 adapter: ${ADAPTER2}"
    fi
  fi

  # Barcode removal (R1)
  if [[ ${BC1_ENABLED} -eq 1 && ${BC1_LENGTH} -gt 0 ]]; then
    if [[ "${BC1_LOCATION}" == "3" ]]; then
      CUTADAPT_CMD+=("-u" "-${BC1_LENGTH}")
      echo "PREP | CUTADAPT | R1 barcode: ${BC1_LENGTH}bp from 3' end"
    else
      CUTADAPT_CMD+=("-u" "${BC1_LENGTH}")
      echo "PREP | CUTADAPT | R1 barcode: ${BC1_LENGTH}bp from 5' end"
    fi
  fi

  # Barcode removal (R2)
  if [[ "${MODE}" == "PE" && ${BC2_ENABLED} -eq 1 && ${BC2_LENGTH} -gt 0 ]]; then
    if [[ "${BC2_LOCATION}" == "3" ]]; then
      CUTADAPT_CMD+=("-U" "-${BC2_LENGTH}")
      echo "PREP | CUTADAPT | R2 barcode: ${BC2_LENGTH}bp from 3' end"
    else
      CUTADAPT_CMD+=("-U" "${BC2_LENGTH}")
      echo "PREP | CUTADAPT | R2 barcode: ${BC2_LENGTH}bp from 5' end"
    fi
  fi

  # Run cutadapt
  echo "PREP | CUTADAPT | Processing reads..."
  if [[ "${MODE}" == "PE" ]]; then
    cutadapt "${CUTADAPT_CMD[@]}" \
             -o preumi_R1.fastq \
             -p preumi_R2.fastq \
             "${R1}" "${R2}" | tee cutadapt_report.txt
  else
    cutadapt "${CUTADAPT_CMD[@]}" \
             -o preumi_R1.fastq \
             "${R1}" | tee cutadapt_report.txt
  fi

  echo "PREP | CUTADAPT | Trimming complete"

  ###########################################################################
  # 6) UMI EXTRACTION
  ###########################################################################

  if [[ ${UMI_ENABLED} -eq 1 && ${UMI_LENGTH} -gt 0 ]]; then
    echo "PREP | UMI | Extracting UMI sequences..."
    echo "PREP | UMI | Length: ${UMI_LENGTH}bp"
    echo "PREP | UMI | Location: ${UMI_LOCATION}' end"

    # Build UMI pattern (e.g., NNNNNN for 6bp UMI)
    UMI_PATTERN=$(printf 'N%.0s' $(seq 1 ${UMI_LENGTH}))
    
    # Location flag
    UMI_END_FLAG=""
    if [[ "${UMI_LOCATION}" == "3" ]]; then
      UMI_END_FLAG="--3prime"
    fi

    # Extract UMI
    if [[ "${MODE}" == "PE" ]]; then
      umi_tools extract \
        --bc-pattern="${UMI_PATTERN}" ${UMI_END_FLAG} \
        -I preumi_R1.fastq -S final_R1.fastq \
        --read2-in preumi_R2.fastq --read2-out final_R2.fastq \
        --log=umi_extract.log
    else
      umi_tools extract \
        --bc-pattern="${UMI_PATTERN}" ${UMI_END_FLAG} \
        -I preumi_R1.fastq -S final_R1.fastq \
        --log=umi_extract.log
      
      # Create empty R2 for SE (consistent tuple shape)
      : > final_R2.fastq
    fi

    echo "PREP | UMI | Extraction complete"

  else
    echo "PREP | UMI | UMI extraction disabled, passing through..."
    
    if [[ "${MODE}" == "PE" ]]; then
      cp -f preumi_R1.fastq final_R1.fastq
      cp -f preumi_R2.fastq final_R2.fastq
    else
      cp -f preumi_R1.fastq final_R1.fastq
      : > final_R2.fastq  # Empty stub for SE
    fi

    echo "PREP | UMI | Pass-through complete"
  fi

  ###########################################################################
  # 7) QC ON FINAL READS
  ###########################################################################

  mkdir -p fastqc_final

  if [[ ${QC_ENABLED} -eq 1 ]]; then
    echo "PREP | QC-FINAL | Running FastQC on final reads..."
    
    if [[ "${MODE}" == "PE" ]]; then
      fastqc --quiet --threads "${THREADS}" -o fastqc_final final_R1.fastq final_R2.fastq
    else
      fastqc --quiet --threads "${THREADS}" -o fastqc_final final_R1.fastq
    fi
    
    echo "PREP | QC-FINAL | FastQC reports generated"
  else
    echo "PREP | QC-FINAL | Skipped (QC disabled)"
    echo "QC disabled." > fastqc_final/README.txt
  fi

  ###########################################################################
  # 8) COMPUTE STATISTICS
  ###########################################################################

  echo "PREP | STATS | Computing read statistics..."

  # Count reads at each stage
  R1_IN=$(count_reads "${R1}")
  R1_PRE=$(count_reads preumi_R1.fastq)
  R1_FINAL=$(count_reads final_R1.fastq)

  if [[ "${MODE}" == "PE" ]]; then
    R2_IN=$(count_reads "${R2}")
    R2_PRE=$(count_reads preumi_R2.fastq)
    R2_FINAL=$(count_reads final_R2.fastq)
  else
    R2_IN=0
    R2_PRE=0
    R2_FINAL=0
  fi

  echo "PREP | STATS | R1: ${R1_IN} → ${R1_PRE} → ${R1_FINAL} reads"
  if [[ "${MODE}" == "PE" ]]; then
    echo "PREP | STATS | R2: ${R2_IN} → ${R2_PRE} → ${R2_FINAL} reads"
  fi

  # Trimming statistics
  {
    echo -e "sample_id\tread\treads_in\treads_preumi\tpct_kept_preumi"
    
    PCT_R1=$(awk -v a="${R1_IN}" -v b="${R1_PRE}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
    printf "%s\tR1\t%s\t%s\t%.2f\n" "${SAMPLE_ID}" "${R1_IN}" "${R1_PRE}" "${PCT_R1}"
    
    if [[ "${MODE}" == "PE" ]]; then
      PCT_R2=$(awk -v a="${R2_IN}" -v b="${R2_PRE}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
      printf "%s\tR2\t%s\t%s\t%.2f\n" "${SAMPLE_ID}" "${R2_IN}" "${R2_PRE}" "${PCT_R2}"
    fi
  } > trim_stats.tsv

  echo "PREP | STATS | Trimming statistics saved: trim_stats.tsv"

  # UMI statistics (if UMI enabled)
  if [[ ${UMI_ENABLED} -eq 1 && ${UMI_LENGTH} -gt 0 ]]; then
    {
      echo -e "sample_id\tumi_enabled\tumi_length\tumi_location\tread\treads_preumi\treads_final\tpct_kept_final"
      
      PCT_R1_UMI=$(awk -v a="${R1_PRE}" -v b="${R1_FINAL}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
      printf "%s\t1\t%s\t%s\tR1\t%s\t%s\t%.2f\n" "${SAMPLE_ID}" "${UMI_LENGTH}" "${UMI_LOCATION}" "${R1_PRE}" "${R1_FINAL}" "${PCT_R1_UMI}"
      
      if [[ "${MODE}" == "PE" ]]; then
        PCT_R2_UMI=$(awk -v a="${R2_PRE}" -v b="${R2_FINAL}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
        printf "%s\t1\t%s\t%s\tR2\t%s\t%s\t%.2f\n" "${SAMPLE_ID}" "${UMI_LENGTH}" "${UMI_LOCATION}" "${R2_PRE}" "${R2_FINAL}" "${PCT_R2_UMI}"
      fi
    } > umi_stats.tsv
    
    echo "PREP | STATS | UMI statistics saved: umi_stats.tsv"
  fi

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "PREP | README | Creating documentation..."

  cat > README_01_trimmed_fastq.txt <<'DOCEOF'
================================================================================
FASTQ PREPROCESSING — !{sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  This directory contains preprocessed FASTQ files ready for alignment.
  
  Processing Pipeline:
    Raw FASTQ → Cutadapt (adapters + barcodes + length) → UMI extraction → Final

FILES
────────────────────────────────────────────────────────────────────────────
  final_R1.fastq              — Cleaned R1 reads (ready for alignment)
  final_R2.fastq              — Cleaned R2 reads (or empty stub for SE)
  
  fastqc_raw/                 — Quality control reports on raw reads
  fastqc_final/               — Quality control reports on final reads
  
  trim_stats.tsv              — Read counts through trimming stages
  umi_stats.tsv               — UMI extraction statistics (if enabled)
  
  cutadapt_report.txt         — Detailed cutadapt trimming report
  umi_extract.log             — UMI extraction log (if enabled)
  preprocess_reads.log        — Complete processing log

PARAMETERS
────────────────────────────────────────────────────────────────────────────
  Sample ID:                  !{sample_id}
  Mode:                       !{(data_type ?: "SE").toString()}
  Threads:                    !{task.cpus}
  
  Adapter Trimming:           !{params.adapter_trimming?.enabled == true ? "enabled" : "disabled"}
  Barcode Removal (R1):       !{params.barcode?.enabled == true ? "enabled" : "disabled"}
  UMI Extraction:             !{params.umi?.enabled == true ? "enabled" : "disabled"}
  QC Tool:                    FastQC
  QC Enabled:                 !{(params.qc?.enabled == null ? true : params.qc?.enabled) ? "yes" : "no"}

QUALITY METRICS
────────────────────────────────────────────────────────────────────────────
  See trim_stats.tsv for detailed read retention rates through each step.
  QC reports (HTML/JSON) are available in fastqc_raw/ and fastqc_final/.

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  Use final_R1.fastq (and final_R2.fastq for PE) for alignment.
  These files have been:
    ✓ Adapter trimmed (if enabled)
    ✓ Barcode removed (if enabled)
    ✓ UMI extracted and appended to read names (if enabled)
    ✓ Length filtered (minimum insert length enforced)
    ✓ Quality checked

NOTES
────────────────────────────────────────────────────────────────────────────
  • Single-pass processing: Faster than traditional multi-step approaches
  • UMI information (if extracted) is embedded in FASTQ read names
  • Empty final_R2.fastq for SE samples maintains consistent file structure

================================================================================
DOCEOF

  echo "PREP | README | Documentation created"

  ###########################################################################
  # 10) SUMMARY
  ###########################################################################

  # Get final file sizes
  FINAL_R1_SIZE=$(stat -c%s final_R1.fastq 2>/dev/null || stat -f%z final_R1.fastq 2>/dev/null || echo "unknown")
  FINAL_R2_SIZE=$(stat -c%s final_R2.fastq 2>/dev/null || stat -f%z final_R2.fastq 2>/dev/null || echo "unknown")

  echo "────────────────────────────────────────────────────────────────────────"
  echo "PREP | SUMMARY | Processing complete for ${SAMPLE_ID}"
  echo "PREP | SUMMARY | Mode: ${MODE}"
  echo "PREP | SUMMARY | Input reads (R1): ${R1_IN}"
  if [[ "${MODE}" == "PE" ]]; then
    echo "PREP | SUMMARY | Input reads (R2): ${R2_IN}"
  fi
  echo "PREP | SUMMARY | Final reads (R1): ${R1_FINAL} ($(awk -v a="${R1_IN}" -v b="${R1_FINAL}" 'BEGIN{printf "%.1f", (a>0)?(b*100.0/a):0}')%)"
  if [[ "${MODE}" == "PE" ]]; then
    echo "PREP | SUMMARY | Final reads (R2): ${R2_FINAL} ($(awk -v a="${R2_IN}" -v b="${R2_FINAL}" 'BEGIN{printf "%.1f", (a>0)?(b*100.0/a):0}')%)"
  fi
  echo "PREP | SUMMARY | Final R1 size: ${FINAL_R1_SIZE} bytes"
  echo "PREP | SUMMARY | Final R2 size: ${FINAL_R2_SIZE} bytes"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "PREP | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}