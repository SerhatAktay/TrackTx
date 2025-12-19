// ============================================================================
// download_srr.nf — SRA to FASTQ Conversion
// ============================================================================
//
// Purpose:
//   Downloads and converts SRA accessions to FASTQ files
//
// Features:
//   • Intelligent caching - reuses existing FASTQs from previous runs
//   • Multi-threaded conversion with fasterq-dump
//   • Optional compression with pigz
//   • Robust validation with header checks
//   • Checksum generation (MD5 or SHA256)
//   • Per-sample provenance documentation
//
// Workflow:
//   1. Check for cached FASTQs in output directory
//   2. Prefetch SRA data (optional, best-effort)
//   3. Convert with fasterq-dump (multi-threaded)
//   4. Normalize filenames (_1/_2 → _R1/_R2)
//   5. Optional compression with pigz
//   6. Validate FASTQ format
//   7. Generate checksums
//   8. Create README with provenance
//
// Inputs:
//   tuple(sample_id, sra_id, condition, timepoint, replicate, paired_end)
//
// Outputs:
//   ${params.output_dir}/01_trimmed_fastq/
//     ├── <SRR>_R1.fastq[.gz]        — Read 1 FASTQ
//     ├── <SRR>_R2.fastq[.gz]        — Read 2 FASTQ (PE only)
//     ├── <SRR>_R1.fastq.md5         — Checksums
//     ├── <SRR>_R2.fastq.md5         — Checksums (PE only)
//     └── README_fastq.txt           — Provenance documentation
//
// Parameters:
//   params.fastq_gzip       : Compress FASTQs with pigz (default: false)
//   params.sra_tmp          : Temp directory for fasterq-dump
//   params.sra_max_size     : Max prefetch size (default: 200G)
//   params.conda_sra        : Conda environment override
//
// ============================================================================

nextflow.enable.dsl = 2

process download_srr {

  tag        { sra_id }
  label      'conda'
  cache      'deep'
  conda      (params.conda_sra ?: "${projectDir}/envs/tracktx.yaml")
  
  // NOTE: Raw FASTQ files are NOT published here to save disk space (~100GB+)
  // They are intermediate files that get processed by prepare_input
  // Only the final trimmed/processed FASTQs are published by prepare_input
  // Raw files remain in work/ directory for Nextflow caching with -resume
  
  publishDir "${params.output_dir}/01_trimmed_fastq",
             mode: 'copy',
             overwrite: true,
             saveAs: { filename ->
               // Only publish checksums and README, NOT the raw FASTQ files
               // This saves ~100GB+ of redundant storage
               if (filename.endsWith('.md5') || 
                   filename.endsWith('.sha256') || 
                   filename == 'README_fastq.txt') {
                 return filename
               }
               return null  // Don't publish raw FASTQs
             }

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), val(sra_id),
          val(condition), val(timepoint), val(replicate),
          val(paired_end)

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("${sra_id}_R1.fastq{,.gz}"),
          path("${sra_id}_R2.fastq{,.gz}"),
          val(condition), val(timepoint), val(replicate)

    path "${sra_id}*.md5",    optional: true, emit: checksums_md5
    path "${sra_id}*.sha256", optional: true, emit: checksums_sha
    path "README_fastq.txt",                  emit: readme

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

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
  SRA_TMP="!{params.sra_tmp ?: ''}"
  SRA_MAX_SIZE="!{params.sra_max_size ?: '200G'}"
  
  # Output directory for caching
  CACHE_DIR="!{params.output_dir}/01_trimmed_fastq"
  mkdir -p "${CACHE_DIR}"

  echo "SRR | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "SRR | CONFIG | Accession: ${SRR}"
  echo "SRR | CONFIG | Layout: $([ "${IS_PE}" == "true" ] && echo "Paired-end" || echo "Single-end")"
  echo "SRR | CONFIG | Threads: ${THREADS}"
  echo "SRR | CONFIG | Compression: ${COMPRESS_FQ}"
  echo "SRR | CONFIG | Cache directory: ${CACHE_DIR}"

  ###########################################################################
  # 2) VALIDATE TOOLS
  ###########################################################################

  echo "SRR | VALIDATE | Checking required tools..."

  if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "SRR | ERROR | fasterq-dump not found in PATH"
    exit 1
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
  # 3) CHECK FOR CACHED FILES
  ###########################################################################

  echo "SRR | CACHE | Checking for existing FASTQs..."

  CACHE_FOUND=0
  for CACHED in "${CACHE_DIR}/${SRR}_R1.fastq" "${CACHE_DIR}/${SRR}_R1.fastq.gz" \
                "${CACHE_DIR}/${SRR}_R2.fastq" "${CACHE_DIR}/${SRR}_R2.fastq.gz"; do
    if [[ -e "${CACHED}" ]]; then
      FILE_SIZE=$(stat -c%s "${CACHED}" 2>/dev/null || stat -f%z "${CACHED}" 2>/dev/null || echo "unknown")
      echo "SRR | CACHE | Found: $(basename "${CACHED}") (${FILE_SIZE} bytes)"
      ln -sf "${CACHED}" .
      CACHE_FOUND=1
    fi
  done

  if [[ ${CACHE_FOUND} -eq 1 ]]; then
    echo "SRR | CACHE | Reusing cached FASTQs from previous run"
  else
    echo "SRR | CACHE | No cached FASTQs found, will download"
  fi

  ###########################################################################
  # 4) DOWNLOAD SRA DATA
  ###########################################################################

  # Only download if R1 doesn't exist yet
  if [[ ! -e "${SRR}_R1.fastq" && ! -e "${SRR}_R1.fastq.gz" ]]; then
    
    # ── Step 4a: Prefetch (optional, best-effort) ──
    if [[ ${HAVE_PREFETCH} -eq 1 ]]; then
      echo "SRR | PREFETCH | Downloading SRA file for ${SRR}..."
      echo "SRR | PREFETCH | Max size: ${SRA_MAX_SIZE}"
      
      if prefetch -O . --verify yes --max-size "${SRA_MAX_SIZE}" "${SRR}" 2>&1 | \
         grep -v "^|" | grep -v "^2025" | grep -v "^$"; then
        echo "SRR | PREFETCH | Download successful"
      else
        echo "SRR | PREFETCH | Prefetch failed or incomplete, will stream directly"
      fi
    fi

    # ── Step 4b: Convert to FASTQ with fasterq-dump ──
    echo "SRR | CONVERT | Converting SRA to FASTQ format..."
    echo "SRR | CONVERT | Using fasterq-dump with ${THREADS} threads"

    # Setup temp directory if specified
    if [[ -n "${SRA_TMP}" ]]; then
      mkdir -p "${SRA_TMP}"
      echo "SRR | CONVERT | Temp directory: ${SRA_TMP}"
      TEMP_ARG="--temp ${SRA_TMP}"
    else
      echo "SRR | CONVERT | Using default temp directory"
      TEMP_ARG=""
    fi

    # Run fasterq-dump with proper error capture
    echo "SRR | CONVERT | Running fasterq-dump --split-files..."
    
    FASTQ_OUTPUT=$(mktemp)
    FASTQ_EXIT=0
    
    if [[ -n "${TEMP_ARG}" ]]; then
      fasterq-dump --split-files -e "${THREADS}" ${TEMP_ARG} -O . "${SRR}" 2>&1 | tee "${FASTQ_OUTPUT}" || FASTQ_EXIT=$?
    else
      fasterq-dump --split-files -e "${THREADS}" -O . "${SRR}" 2>&1 | tee "${FASTQ_OUTPUT}" || FASTQ_EXIT=$?
    fi
    
    # Check for errors
    if [[ ${FASTQ_EXIT} -ne 0 ]]; then
      echo "SRR | ERROR | fasterq-dump failed with exit code ${FASTQ_EXIT}"
      echo "SRR | ERROR | Full output:"
      cat "${FASTQ_OUTPUT}"
      rm -f "${FASTQ_OUTPUT}"
      exit 1
    fi
    
    # Show summary stats only
    grep -E "spots read|reads read|reads written" "${FASTQ_OUTPUT}" || true
    rm -f "${FASTQ_OUTPUT}"

    echo "SRR | CONVERT | Conversion complete"

    # ── Step 4c: Normalize filenames ──
    echo "SRR | CONVERT | Normalizing FASTQ filenames..."
    
    # Handle different output patterns from fasterq-dump
    if [[ -f "${SRR}_1.fastq" ]]; then
      mv -f "${SRR}_1.fastq" "${SRR}_R1.fastq"
      echo "SRR | CONVERT | Renamed ${SRR}_1.fastq → ${SRR}_R1.fastq"
    fi
    
    if [[ -f "${SRR}_2.fastq" ]]; then
      mv -f "${SRR}_2.fastq" "${SRR}_R2.fastq"
      echo "SRR | CONVERT | Renamed ${SRR}_2.fastq → ${SRR}_R2.fastq"
    fi
    
    if [[ -f "${SRR}.fastq" ]]; then
      mv -f "${SRR}.fastq" "${SRR}_R1.fastq"
      echo "SRR | CONVERT | Renamed ${SRR}.fastq → ${SRR}_R1.fastq"
    fi

  else
    echo "SRR | CONVERT | Skipping download, using cached files"
  fi

  ###########################################################################
  # 5) OPTIONAL COMPRESSION
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

  ###########################################################################
  # 6) VALIDATE OUTPUT FILES
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
    echo "SRR | ERROR | R1 FASTQ file missing or empty"
    exit 1
  fi

  R1_SIZE=$(stat -c%s "${R1_FILE}" 2>/dev/null || stat -f%z "${R1_FILE}" 2>/dev/null || echo "unknown")
  echo "SRR | VALIDATE | R1 file: ${R1_FILE} (${R1_SIZE} bytes)"

  # Quick header check for uncompressed files
  if [[ "${R1_FILE}" == *.fastq ]]; then
    FIRST_LINE=$(head -n 1 "${R1_FILE}" 2>/dev/null || echo "")
    if [[ ! "${FIRST_LINE}" =~ ^@ ]]; then
      echo "SRR | ERROR | R1 FASTQ header validation failed"
      echo "SRR | ERROR | First line does not start with '@': ${FIRST_LINE}"
      exit 1
    fi
    echo "SRR | VALIDATE | R1 FASTQ header looks valid"
  fi

  # Validate R2 (if paired-end)
  if [[ "${IS_PE}" == "true" ]]; then
    if [[ -z "${R2_FILE}" || ! -s "${R2_FILE}" ]]; then
      echo "SRR | ERROR | Paired-end mode but R2 FASTQ missing or empty"
      exit 1
    fi

    R2_SIZE=$(stat -c%s "${R2_FILE}" 2>/dev/null || stat -f%z "${R2_FILE}" 2>/dev/null || echo "unknown")
    echo "SRR | VALIDATE | R2 file: ${R2_FILE} (${R2_SIZE} bytes)"

    # Quick header check for uncompressed R2
    if [[ "${R2_FILE}" == *.fastq ]]; then
      FIRST_LINE=$(head -n 1 "${R2_FILE}" 2>/dev/null || echo "")
      if [[ ! "${FIRST_LINE}" =~ ^@ ]]; then
        echo "SRR | ERROR | R2 FASTQ header validation failed"
        echo "SRR | ERROR | First line does not start with '@': ${FIRST_LINE}"
        exit 1
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
  # 7) GENERATE CHECKSUMS
  ###########################################################################

  echo "SRR | CHECKSUM | Generating checksums..."

  # Generate checksum for R1
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "${R1_FILE}" > "${R1_FILE}.md5"
    echo "SRR | CHECKSUM | ${R1_FILE}.md5 created"
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${R1_FILE}" > "${R1_FILE}.sha256"
    echo "SRR | CHECKSUM | ${R1_FILE}.sha256 created"
  else
    echo "SRR | CHECKSUM | WARNING: No checksum tool available (md5sum/shasum)"
  fi

  # Generate checksum for R2 if exists
  if [[ -n "${R2_FILE}" && -f "${R2_FILE}" ]]; then
    if command -v md5sum >/dev/null 2>&1; then
      md5sum "${R2_FILE}" > "${R2_FILE}.md5"
      echo "SRR | CHECKSUM | ${R2_FILE}.md5 created"
    elif command -v shasum >/dev/null 2>&1; then
      shasum -a 256 "${R2_FILE}" > "${R2_FILE}.sha256"
      echo "SRR | CHECKSUM | ${R2_FILE}.sha256 created"
    fi
  fi

  ###########################################################################
  # 8) CREATE README
  ###########################################################################

  echo "SRR | README | Creating provenance documentation..."

  cat > README_fastq.txt <<'DOCEOF'
================================================================================
SRA FASTQ DOWNLOAD — !{sra_id}
================================================================================

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample ID:    !{sample_id}
  SRA Accession: !{sra_id}
  Condition:    !{condition}
  Timepoint:    !{timepoint}
  Replicate:    !{replicate}
  Layout:       !{paired_end == "true" ? "Paired-end" : "Single-end"}
  Compression:  !{(params.fastq_gzip == null) ? "false" : (params.fastq_gzip as boolean ? "true" : "false")}
  Downloaded:   $(date -u +"%Y-%m-%d %H:%M:%S UTC")

PIPELINE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Pipeline:     TrackTx PRO-seq Analysis
  Module:       02_download_srr
  Threads:      !{task.cpus}

FILES
────────────────────────────────────────────────────────────────────────────
  R1: !{sra_id}_R1.fastq!{(params.fastq_gzip == null) ? "" : (params.fastq_gzip as boolean ? ".gz" : "")}
  R2: !{paired_end == "true" ? sra_id + "_R2.fastq" + ((params.fastq_gzip == null) ? "" : (params.fastq_gzip as boolean ? ".gz" : "")) : "N/A (single-end)"}

CACHING
────────────────────────────────────────────────────────────────────────────
  These files are cached in the output directory for reuse across pipeline runs.
  Nextflow will automatically detect and reuse these files with -resume.

NOTES
────────────────────────────────────────────────────────────────────────────
  • Files downloaded using SRA Toolkit (fasterq-dump)
  • Checksums (MD5 or SHA256) generated for data integrity
  • For access-controlled data (dbGaP), ensure SRA Toolkit is properly configured:
      vdb-config --interactive
      vdb-config --import <project.krt>

================================================================================
DOCEOF

  echo "SRR | README | Documentation created"

  ###########################################################################
  # 9) SUMMARY
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