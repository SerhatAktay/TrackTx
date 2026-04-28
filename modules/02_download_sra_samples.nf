// ============================================================================
// download_sra_samples.nf — SRA to FASTQ Conversion
// ============================================================================
//
// Purpose:
//   Downloads and converts SRA accessions to FASTQ files
//
// Features:
//   • storeDir caching — outputs persist in results/00_sra_cache/<SRR>/ across work/ deletions
//   • Multi-threaded conversion with fasterq-dump
//   • Optional compression with pigz
//   • Robust validation with header checks
//
// Workflow:
//   1. Prefetch SRA data (optional, best-effort)
//   2. Convert with fasterq-dump (multi-threaded)
//   3. Normalize filenames (_1/_2 → _R1/_R2)
//   4. Optional compression with pigz
//   5. Validate FASTQ format
//
// Inputs:
//   tuple(sample_id, sra_id, condition, timepoint, replicate, paired_end)
//
// Outputs:
//   When publish_sra_fastq=true (default):
//     ${params.output_dir}/00_sra_cache/<SRR>/
//       ├── <SRR>_R1.fastq[.gz]      — Read 1 FASTQ  (persisted as storeDir cache)
//       └── <SRR>_R2.fastq[.gz]      — Read 2 FASTQ (PE only)
//   When publish_sra_fastq=false:
//     ${params.output_dir}/.sra_cache/<SRR>/  (hidden; saves ~100 GB+ in results/)
//
// Parameters:
//   params.publish_sra_fastq : Save raw FASTQs to results/00_sra_cache/ (default: true)
//                              false → cache to results/.sra_cache/ (hidden from results view)
//   params.fastq_gzip        : Compress FASTQs with pigz (default: false)
//   params.sra_tmp           : Temp directory for fasterq-dump
//   params.sra_max_size      : Max prefetch size (default: 200G)
//   params.conda_sra         : Conda environment override
//
// ============================================================================

nextflow.enable.dsl = 2

process download_sra_samples {

  tag        { sra_id }
  label      'conda'
  cache      'deep'
  conda      (params.conda_sra ?: "${projectDir}/envs/tracktx.yaml")
  
  // storeDir persists raw FASTQs so that re-runs skip the download entirely — even after
  // the work/ directory has been deleted. Cache is keyed per SRA accession so adding new
  // samples never invalidates existing downloads.
  //
  // publish_sra_fastq=true  (default): cache in results/00_sra_cache/<SRR>/  — visible in results
  // publish_sra_fastq=false          : cache in results/.sra_cache/<SRR>/    — hidden from results
  //   In false mode the primary re-run safeguard is the trimmed-FASTQ check in main.nf (STEP 3).
  //   The hidden .sra_cache/ acts as a secondary fallback. Clean it up with: rm -rf results/.sra_cache/
  storeDir { params.get('publish_sra_fastq') != false
      ? "${params.output_dir}/00_sra_cache/${sra_id}"
      : "${params.output_dir}/.sra_cache/${sra_id}" }

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), val(sra_id),
          val(condition), val(timepoint), val(replicate),
          val(paired_end)

  // ── Outputs ───────────────────────────────────────────────────────────────
  // We always emit _R1.fastq and _R2.fastq. When compression is on, the script
  // creates symlinks .fastq -> .fastq.gz so downstream receives consistent paths.
  // NOTE: preprocess_and_quality_filter_reads detects gzipped content with .fastq extension and renames
  // to .gz so FastQC/cutadapt get correct format (they detect by extension).
  output:
    tuple val(sample_id),
          path("${sra_id}_R1.fastq"),
          path("${sra_id}_R2.fastq"),
          val(condition), val(timepoint), val(replicate)

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  sraTmp        = (params.sra_tmp ?: '').toString()
  sraMaxSize    = (params.sra_max_size ?: '200G').toString()
  sraSource     = (params.sra_download_source ?: 'auto').toString()
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  tracktx_error() {
    local module="\$1" problem="\$2" fix="\$3" code="\${4:-1}"
    echo "" >&2
    echo "═══════════════════════════════════════════════════════════════════════" >&2
    echo "TRACKTX ERROR" >&2
    echo "═══════════════════════════════════════════════════════════════════════" >&2
    echo "Module:  \${module}" >&2
    echo "Problem: \${problem}" >&2
    echo "Fix:     \${fix}" >&2
    echo "═══════════════════════════════════════════════════════════════════════" >&2
    exit "\$code"
  }

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "SRR | START | sample=!{sample_id} | accession=!{sra_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sample_id}"
  SRR="!{sra_id}"
  CONDITION="!{condition}"
  TIMEPOINT="!{timepoint}"
  REPLICATE="!{replicate}"
  IS_PE="!{paired_end}"
  THREADS=!{task.cpus}
  
  # Compression settings
  COMPRESS_FQ="!{(params.fastq_gzip == null) ? 'false' : (params.fastq_gzip as boolean ? 'true' : 'false')}"
  
  # SRA settings
  SRA_TMP="!{sraTmp}"
  SRA_MAX_SIZE="!{sraMaxSize}"
  SRA_SOURCE="!{sraSource}"
  
  echo "SRR | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "SRR | CONFIG | Accession: ${SRR}"
  echo "SRR | CONFIG | Layout: $([ "${IS_PE}" == "true" ] && echo "Paired-end" || echo "Single-end")"
  echo "SRR | CONFIG | Threads: ${THREADS}"
  echo "SRR | CONFIG | Compression: ${COMPRESS_FQ}"

  ###########################################################################
  # 2) VALIDATE TOOLS
  ###########################################################################

  echo "SRR | VALIDATE | Checking required tools..."

  if ! command -v fasterq-dump >/dev/null 2>&1; then
    tracktx_error "download_sra_samples" "fasterq-dump not found in PATH" "Install sra-tools or use -profile docker"
  fi
  echo "SRR | VALIDATE | fasterq-dump: $(which fasterq-dump)"

  HAVE_PREFETCH=0
  if command -v prefetch >/dev/null 2>&1; then
    HAVE_PREFETCH=1
    echo "SRR | VALIDATE | prefetch: $(which prefetch)"
  else
    echo "SRR | VALIDATE | prefetch not available (will stream directly)"
  fi

  if command -v pigz >/dev/null 2>&1; then
    echo "SRR | VALIDATE | pigz: $(which pigz) (for compression)"
  else
    echo "SRR | VALIDATE | pigz not available (compression disabled)"
  fi

  ###########################################################################
  # 3) DOWNLOAD SRA DATA
  ###########################################################################

  # Only download if R1 doesn't exist yet
  if [[ ! -e "${SRR}_R1.fastq" && ! -e "${SRR}_R1.fastq.gz" ]]; then
    
    TRY_ENA=0
    if [[ "${SRA_SOURCE}" == "ena" ]]; then
      TRY_ENA=1
      echo "SRR | SOURCE | Using ENA (European Nucleotide Archive) for download"
    fi
    
    # ── Try NCBI (unless ena-only) ──
    NCBI_OK=0
    if [[ ${TRY_ENA} -eq 0 ]]; then
      if [[ ${HAVE_PREFETCH} -eq 1 ]]; then
        echo "SRR | PREFETCH | Downloading SRA file for ${SRR}..."
        prefetch -O . --verify yes --max-size "${SRA_MAX_SIZE}" "${SRR}" 2>&1 | grep -v "^|" | grep -v "^202" | grep -v "^$" || true
      fi
      echo "SRR | CONVERT | Converting SRA to FASTQ (fasterq-dump)..."
      FASTQ_OUTPUT=\$(mktemp)
      TEMP_ARG=""
      [[ -n "${SRA_TMP}" ]] && mkdir -p "${SRA_TMP}" && TEMP_ARG="--temp ${SRA_TMP}"
      if [[ -n "\${TEMP_ARG}" ]]; then
        fasterq-dump --split-files -e "${THREADS}" \${TEMP_ARG} -O . "${SRR}" 2>&1 | tee "\${FASTQ_OUTPUT}" || true
      else
        fasterq-dump --split-files -e "${THREADS}" -O . "${SRR}" 2>&1 | tee "\${FASTQ_OUTPUT}" || true
      fi
      grep -E "spots read|reads read|reads written" "\${FASTQ_OUTPUT}" 2>/dev/null || true
      rm -f "\${FASTQ_OUTPUT}"
      [[ -f "${SRR}_1.fastq" || -f "${SRR}_2.fastq" || -f "${SRR}.fastq" ]] && NCBI_OK=1
    fi
    
    # ── Fallback to ENA if NCBI failed (or ena-only) ──
    if [[ ${NCBI_OK} -eq 0 && ( ${TRY_ENA} -eq 1 || "${SRA_SOURCE}" == "auto" ) ]]; then
      echo "SRR | ENA | Downloading from ENA (NCBI unavailable or ena-only mode)..."
      ENA_REPORT=\$(mktemp)
      if curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp" -o "\${ENA_REPORT}" 2>/dev/null; then
        ENA_FTP=\$(tail -n 1 "\${ENA_REPORT}" | cut -f2)
        rm -f "\${ENA_REPORT}"
        if [[ -n "\${ENA_FTP}" && "\${ENA_FTP}" != "fastq_ftp" ]]; then
          # Launch all file downloads in parallel (R1 + R2 simultaneously)
          ENA_PIDS=()
          IX=0
          for URL in \$(echo "\${ENA_FTP}" | tr ';' ' '); do
            [[ -z "\${URL}" ]] && continue
            IX=$((IX+1))
            FNAME="\${URL##*/}"
            FURL="https://\${URL}"
            echo "SRR | ENA | Starting download \${IX}: \${FNAME}"
            {
              if curl -sfL "\${FURL}" -o "\${FNAME}" 2>/dev/null; then
                echo "SRR | ENA | Downloaded \${FNAME}"
              else
                echo "SRR | ENA | HTTPS failed for \${FNAME}, trying FTP..."
                curl -sfL "ftp://\${URL}" -o "\${FNAME}"
              fi
            } &
            ENA_PIDS+=(\$!)
          done
          ENA_FAILED=0
          for pid in "\${ENA_PIDS[@]}"; do wait "\${pid}" || ENA_FAILED=1; done
          if [[ \${ENA_FAILED} -ne 0 ]]; then
            echo "SRR | ENA | One or more downloads failed"
          fi
          NCBI_OK=1
        else
          echo "SRR | ENA | No FASTQ URLs in ENA report"
        fi
      else
        rm -f "\${ENA_REPORT}"
        echo "SRR | ENA | Failed to fetch ENA filereport"
      fi
    fi
    
    # ── Normalize filenames ──
    echo "SRR | CONVERT | Normalizing FASTQ filenames..."
    if [[ -f "${SRR}_1.fastq" ]]; then
      mv -f "${SRR}_1.fastq" "${SRR}_R1.fastq"
    elif [[ -f "${SRR}_1.fastq.gz" ]]; then
      mv -f "${SRR}_1.fastq.gz" "${SRR}_R1.fastq.gz"
    fi
    if [[ -f "${SRR}_2.fastq" ]]; then
      mv -f "${SRR}_2.fastq" "${SRR}_R2.fastq"
    elif [[ -f "${SRR}_2.fastq.gz" ]]; then
      mv -f "${SRR}_2.fastq.gz" "${SRR}_R2.fastq.gz"
    fi
    if [[ -f "${SRR}.fastq" ]]; then
      mv -f "${SRR}.fastq" "${SRR}_R1.fastq"
    elif [[ -f "${SRR}.fastq.gz" ]]; then
      mv -f "${SRR}.fastq.gz" "${SRR}_R1.fastq.gz"
    fi
    
    # Validate we got files
    if [[ ! -e "${SRR}_R1.fastq" && ! -e "${SRR}_R1.fastq.gz" ]]; then
      tracktx_error "download_sra_samples" "Download failed from both NCBI and ENA" "Check SRA accession and network"
    fi
    echo "SRR | CONVERT | Conversion complete"

  else
    echo "SRR | CONVERT | Skipping download, using cached files"
  fi

  ###########################################################################
  # 4) OPTIONAL COMPRESSION
  ###########################################################################

  if [[ "${COMPRESS_FQ}" == "true" ]]; then
    echo "SRR | COMPRESS | Compressing FASTQs with pigz..."
    
    if command -v pigz >/dev/null 2>&1; then
      for FASTQ in "${SRR}_R1.fastq" "${SRR}_R2.fastq"; do
        if [[ -f "${FASTQ}" ]]; then
          echo "SRR | COMPRESS | Compressing ${FASTQ}..."
          pigz -p "${THREADS}" --force "${FASTQ}"
          echo "SRR | COMPRESS | Created ${FASTQ}.gz"
        fi
      done
    else
      echo "SRR | COMPRESS | WARNING: pigz not found, leaving FASTQs uncompressed"
    fi
  else
    echo "SRR | COMPRESS | Compression disabled, FASTQs will remain uncompressed"
  fi

  # Ensure .fastq outputs exist for Nextflow (when compressed, symlink .fastq -> .fastq.gz)
  for F in "${SRR}_R1.fastq" "${SRR}_R2.fastq"; do
    if [[ -f "${F}.gz" && ! -e "${F}" ]]; then
      ln -sf "$(basename "${F}.gz")" "${F}"
      echo "SRR | OUTPUT | Created symlink ${F} -> ${F}.gz"
    fi
  done

  ###########################################################################
  # 5) VALIDATE OUTPUT FILES
  ###########################################################################

  echo "SRR | VALIDATE | Validating output FASTQ files..."

  # Determine actual filenames based on compression
  R1_FILE=""
  R2_FILE=""

  if [[ -f "${SRR}_R1.fastq.gz" ]]; then
    R1_FILE="${SRR}_R1.fastq.gz"
  elif [[ -f "${SRR}_R1.fastq" ]]; then
    R1_FILE="${SRR}_R1.fastq"
  fi

  if [[ -f "${SRR}_R2.fastq.gz" ]]; then
    R2_FILE="${SRR}_R2.fastq.gz"
  elif [[ -f "${SRR}_R2.fastq" ]]; then
    R2_FILE="${SRR}_R2.fastq"
  fi

  # Validate R1 (required)
  if [[ -z "${R1_FILE}" || ! -s "${R1_FILE}" ]]; then
    tracktx_error "download_sra_samples" "R1 FASTQ file missing or empty" "Check fasterq-dump output"
  fi

  R1_SIZE=$(stat -c%s "${R1_FILE}" 2>/dev/null || stat -f%z "${R1_FILE}" 2>/dev/null || echo "unknown")
  echo "SRR | VALIDATE | R1 file: ${R1_FILE} (${R1_SIZE} bytes)"

  # Quick header check for uncompressed files
  if [[ "${R1_FILE}" == *.fastq ]]; then
    FIRST_LINE=$(head -n 1 "${R1_FILE}" 2>/dev/null || echo "")
    if [[ ! "${FIRST_LINE}" =~ ^@ ]]; then
      tracktx_error "download_sra_samples" "R1 FASTQ header validation failed (first line: ${FIRST_LINE})" "Check SRA data integrity"
    fi
    echo "SRR | VALIDATE | R1 FASTQ header looks valid"
  fi

  # Validate R2 (if paired-end)
  if [[ "${IS_PE}" == "true" ]]; then
    if [[ -z "${R2_FILE}" || ! -s "${R2_FILE}" ]]; then
      tracktx_error "download_sra_samples" "Paired-end mode but R2 FASTQ missing or empty" "Check SRA layout and fasterq-dump"
    fi

    R2_SIZE=$(stat -c%s "${R2_FILE}" 2>/dev/null || stat -f%z "${R2_FILE}" 2>/dev/null || echo "unknown")
    echo "SRR | VALIDATE | R2 file: ${R2_FILE} (${R2_SIZE} bytes)"

    # Quick header check for uncompressed R2
    if [[ "${R2_FILE}" == *.fastq ]]; then
      FIRST_LINE=$(head -n 1 "${R2_FILE}" 2>/dev/null || echo "")
      if [[ ! "${FIRST_LINE}" =~ ^@ ]]; then
        tracktx_error "download_sra_samples" "R2 FASTQ header validation failed (first line: ${FIRST_LINE})" "Check SRA data integrity"
      fi
      echo "SRR | VALIDATE | R2 FASTQ header looks valid"
    fi
  else
    echo "SRR | VALIDATE | Single-end mode, R2 not expected"
    # Create empty R2 placeholder for Nextflow output consistency
    touch "${SRR}_R2.fastq"
    echo "SRR | VALIDATE | Created empty R2 placeholder for single-end mode"
  fi

  ###########################################################################
  # 6) SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "SRR | SUMMARY | Download complete for ${SRR}"
  echo "SRR | SUMMARY | R1: ${R1_FILE} (${R1_SIZE} bytes)"
  if [[ -n "${R2_FILE}" ]]; then
    echo "SRR | SUMMARY | R2: ${R2_FILE} (${R2_SIZE} bytes)"
  fi
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "SRR | COMPLETE | sample=${SAMPLE_ID} | accession=${SRR} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}