// ============================================================================
// download_gtf.nf — Genome Annotation Download and Processing
// ============================================================================
//
// Purpose:
//   Downloads reference GTF annotation and generates derived gene catalogs
//
// Features:
//   • Fetches GTF from multiple sources (local, custom URL, or UCSC)
//   • Generates gene catalog tables (genes.tsv, TSS.bed, TES.bed)
//   • Caches artifacts for reuse across pipeline runs
//   • Deterministic outputs with coordinate sorting
//
// Sources (in priority order):
//   1. Local file (if exists in genomes/ directory)
//   2. Custom path (if reference_genome='other' with --gtf_path)
//   3. Custom URL (if reference_genome='other' with --gtf_url)
//   4. UCSC (rsync → https fallback)
//
// Inputs:
//   None (driven by params)
//
//   Required Parameters:
//     params.reference_genome : genome ID (hg38, mm39, etc.) or 'other'
//     params.genome_cache     : cache directory for downloads
//     params.assets_dir       : persistent artifact storage
//
//   Custom GTF (when reference_genome='other'):
//     params.gtf_path : local file path
//     params.gtf_url  : remote URL
//
// Outputs:
//   ${params.output_dir}/00_references/${reference_genome}/
//     ├── ${reference_genome}.gtf        — Full GTF annotation
//     ├── ${reference_genome}.genes.tsv  — Gene catalog table
//     ├── ${reference_genome}.tss.bed    — Transcription start sites (sorted)
//     └── ${reference_genome}.tes.bed    — Transcription end sites (sorted)
//
// Caching:
//   Artifacts stored in: ${assets_dir}/annotation/${reference_genome}/
//   This enables fast reuse across pipeline runs
//
// ============================================================================

nextflow.enable.dsl = 2

process download_gtf {

  tag        { params.reference_genome }
  label      'conda'
  cache      'deep'
  conda      'envs/tracktx.yaml'

  // Persistent storage for caching across runs
  storeDir   "${params.assets_dir ?: "${projectDir}/assets"}/annotation/${params.reference_genome}"
  
  // Lightweight links in output directory
  publishDir "${params.output_dir}/00_references/${params.reference_genome}", 
             mode: 'link', 
             overwrite: true

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    path("${params.reference_genome}.gtf")
    path("${params.reference_genome}.genes.tsv"), emit: genes
    path("${params.reference_genome}.tss.bed"),   emit: tss
    path("${params.reference_genome}.tes.bed"),   emit: tes

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "GTF | START | genome=!{params.reference_genome} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  ASM='!{params.reference_genome}'
  CACHE_DIR='!{params.genome_cache ?: "/tmp/genomes_cache"}'
  mkdir -p "${CACHE_DIR}"

  # Output file names
  OUT_GTF="${ASM}.gtf"
  OUT_GENES="${ASM}.genes.tsv"
  OUT_TSS="${ASM}.tss.bed"
  OUT_TES="${ASM}.tes.bed"

  # Cache file paths
  CACHE_GTF="${CACHE_DIR}/${ASM}.gtf"
  CACHE_GENES="${CACHE_DIR}/${ASM}.genes.tsv"
  CACHE_TSS="${CACHE_DIR}/${ASM}.tss.bed"
  CACHE_TES="${CACHE_DIR}/${ASM}.tes.bed"

  # Custom GTF sources
  CUSTOM_PATH='!{params.gtf_path ?: ""}'
  CUSTOM_URL='!{params.gtf_url ?: ""}'

  # UCSC download URLs
  RSYNC_URL="rsync://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/bigZips/genes/${ASM}.ncbiRefSeq.gtf.gz"
  HTTPS_URL="https://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/bigZips/genes/${ASM}.ncbiRefSeq.gtf.gz"

  echo "GTF | CONFIG | Assembly: ${ASM}"
  echo "GTF | CONFIG | Cache directory: ${CACHE_DIR}"
  echo "GTF | CONFIG | Custom path: ${CUSTOM_PATH:-none}"
  echo "GTF | CONFIG | Custom URL: ${CUSTOM_URL:-none}"

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  else
    PYTHON_CMD="python3"
  fi

  # Validate dependencies
  if ! ${PYTHON_CMD} --version >/dev/null 2>&1; then
    echo "GTF | ERROR | Python not found (tried: ${PYTHON_CMD})"
    exit 1
  fi
  echo "GTF | CONFIG | Python: $(${PYTHON_CMD} --version)"

  ###########################################################################
  # 2) FAST PATH - Use Complete Cache
  ###########################################################################

  if [[ -s "${CACHE_GTF}" && -s "${CACHE_GENES}" && \
        -s "${CACHE_TSS}" && -s "${CACHE_TES}" ]]; then
    echo "GTF | CACHE | Found complete cached annotation set"
    echo "GTF | CACHE | Reusing: ${CACHE_GTF}"
    echo "GTF | CACHE | Reusing: ${CACHE_GENES}"
    echo "GTF | CACHE | Reusing: ${CACHE_TSS}"
    echo "GTF | CACHE | Reusing: ${CACHE_TES}"
    
    cp -f "${CACHE_GTF}"   "${OUT_GTF}"
    cp -f "${CACHE_GENES}" "${OUT_GENES}"
    cp -f "${CACHE_TSS}"   "${OUT_TSS}"
    cp -f "${CACHE_TES}"   "${OUT_TES}"
    
    TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    echo "════════════════════════════════════════════════════════════════════════"
    echo "GTF | COMPLETE | Using cached files | ts=${TIMESTAMP_END}"
    echo "════════════════════════════════════════════════════════════════════════"
    exit 0
  fi

  ###########################################################################
  # 3) FETCH GTF (if not cached)
  ###########################################################################

  if [[ -s "${CACHE_GTF}" ]]; then
    echo "GTF | CACHE | GTF found in cache, will regenerate derived files"
    cp -f "${CACHE_GTF}" "${OUT_GTF}"
  else
    echo "GTF | FETCH | GTF not in cache, downloading..."
    
    # Create temporary directory for download
    TEMP_DIR="$(mktemp -d "${CACHE_DIR}/.gtf_${ASM}.XXXXXX")"
    trap 'rm -rf "${TEMP_DIR}"' EXIT
    GTF_TEMP="${TEMP_DIR}/${ASM}.gtf"
    
    FETCH_SUCCESS=0

    # ── Source 1: Local file in genomes/ directory ──
    echo "GTF | FETCH | Checking for local GTF file..."
    LOCAL_GTF_1="${CACHE_DIR}/local_genomes/${ASM}.gtf"
    LOCAL_GTF_2="${CACHE_DIR}/local_genomes/${ASM}/${ASM}.gtf"
    
    if [[ -s "${LOCAL_GTF_1}" ]]; then
      echo "GTF | FETCH | Found local GTF: ${LOCAL_GTF_1}"
      cp -f "${LOCAL_GTF_1}" "${GTF_TEMP}"
      FETCH_SUCCESS=1
    elif [[ -s "${LOCAL_GTF_2}" ]]; then
      echo "GTF | FETCH | Found local GTF: ${LOCAL_GTF_2}"
      cp -f "${LOCAL_GTF_2}" "${GTF_TEMP}"
      FETCH_SUCCESS=1
    fi

    # ── Source 2: Custom path (for reference_genome='other') ──
    if [[ ${FETCH_SUCCESS} -eq 0 && "${ASM}" == "other" && -n "${CUSTOM_PATH}" ]]; then
      if [[ -e "${CUSTOM_PATH}" ]]; then
        echo "GTF | FETCH | Using custom GTF path: ${CUSTOM_PATH}"
        
        if [[ "${CUSTOM_PATH}" == *.gz ]]; then
          echo "GTF | FETCH | Decompressing gzipped GTF..."
          gzip -cd "${CUSTOM_PATH}" > "${GTF_TEMP}"
        else
          cp -f "${CUSTOM_PATH}" "${GTF_TEMP}"
        fi
        FETCH_SUCCESS=1
      else
        echo "GTF | ERROR | Custom GTF path does not exist: ${CUSTOM_PATH}"
        exit 1
      fi
    fi

    # ── Source 3: Custom URL (for reference_genome='other') ──
    if [[ ${FETCH_SUCCESS} -eq 0 && "${ASM}" == "other" && -n "${CUSTOM_URL}" ]]; then
      echo "GTF | FETCH | Downloading from custom URL: ${CUSTOM_URL}"
      
      if [[ "${CUSTOM_URL}" == *.gz ]]; then
        echo "GTF | FETCH | Downloading and decompressing..."
        if curl -fsSL --retry 3 --retry-delay 4 "${CUSTOM_URL}" | gzip -cd > "${GTF_TEMP}"; then
          FETCH_SUCCESS=1
          echo "GTF | FETCH | Download successful"
        else
          echo "GTF | ERROR | Failed to download from: ${CUSTOM_URL}"
        fi
      else
        if curl -fsSL --retry 3 --retry-delay 4 "${CUSTOM_URL}" > "${GTF_TEMP}"; then
          FETCH_SUCCESS=1
          echo "GTF | FETCH | Download successful"
        else
          echo "GTF | ERROR | Failed to download from: ${CUSTOM_URL}"
        fi
      fi
    fi

    # ── Source 4: UCSC (rsync → https fallback) ──
    if [[ ${FETCH_SUCCESS} -eq 0 && "${ASM}" != "other" ]]; then
      echo "GTF | FETCH | Downloading from UCSC: ${ASM}"
      
      # Try rsync first (faster)
      if command -v rsync >/dev/null 2>&1; then
        echo "GTF | FETCH | Attempting rsync download..."
        if rsync --quiet "${RSYNC_URL}" - 2>/dev/null | gunzip -c > "${GTF_TEMP}"; then
          FETCH_SUCCESS=1
          echo "GTF | FETCH | Rsync download successful"
        else
          echo "GTF | FETCH | Rsync failed, will try HTTPS..."
        fi
      else
        echo "GTF | FETCH | Rsync not available, using HTTPS..."
      fi
      
      # Fallback to HTTPS
      if [[ ${FETCH_SUCCESS} -eq 0 ]]; then
        echo "GTF | FETCH | Downloading via HTTPS..."
        if curl -fsSL --retry 3 --retry-delay 4 "${HTTPS_URL}" | gunzip -c > "${GTF_TEMP}"; then
          FETCH_SUCCESS=1
          echo "GTF | FETCH | HTTPS download successful"
        else
          echo "GTF | ERROR | HTTPS download failed"
        fi
      fi
    fi

    # Validate download
    if [[ ${FETCH_SUCCESS} -eq 0 || ! -s "${GTF_TEMP}" ]]; then
      echo "GTF | ERROR | Failed to obtain GTF for assembly: ${ASM}"
      echo "GTF | ERROR | Tried: local file, custom path/URL, UCSC"
      exit 1
    fi

    # Validate GTF format (non-fatal warning)
    echo "GTF | VALIDATE | Checking GTF format..."
    if ! awk 'BEGIN{FS="\\t"} NF>=9{ok=1; exit} END{exit (ok?0:1)}' "${GTF_TEMP}"; then
      echo "GTF | WARNING | GTF may not have standard 9-column format"
      echo "GTF | WARNING | Proceeding anyway..."
    else
      echo "GTF | VALIDATE | GTF format appears valid"
    fi

    # Get file size for reporting
    GTF_SIZE=$(stat -f%z "${GTF_TEMP}" 2>/dev/null || stat -c%s "${GTF_TEMP}" 2>/dev/null || echo "unknown")
    echo "GTF | FETCH | Downloaded GTF size: ${GTF_SIZE} bytes"

    # Move to cache
    mv -f "${GTF_TEMP}" "${CACHE_GTF}"
    cp -f "${CACHE_GTF}" "${OUT_GTF}"
    echo "GTF | FETCH | GTF saved to cache: ${CACHE_GTF}"
  fi

  ###########################################################################
  # 4) GENERATE DERIVED CATALOGS
  ###########################################################################

  echo "GTF | PROCESS | Generating gene catalogs from GTF..."
  
  # Create temporary work directory
  WORK_DIR="$(mktemp -d "${CACHE_DIR}/.build_${ASM}.XXXXXX")"
  trap 'rm -rf "${WORK_DIR}"' EXIT

  # Locate Python script
  SCRIPT_PATH="!{projectDir}/bin/gtf_to_catalog.py"
  if [[ ! -f "${SCRIPT_PATH}" ]]; then
    # Fallback to PATH if not found (for container environments)
    SCRIPT_PATH="gtf_to_catalog.py"
    echo "GTF | PROCESS | Using gtf_to_catalog.py from PATH"
  else
    echo "GTF | PROCESS | Using script: ${SCRIPT_PATH}"
  fi

  # Generate catalogs
  echo "GTF | PROCESS | Running gtf_to_catalog.py..."
  ${PYTHON_CMD} "${SCRIPT_PATH}" \
    "${CACHE_GTF}" \
    "${WORK_DIR}/genes.tsv" \
    "${WORK_DIR}/tss.bed" \
    "${WORK_DIR}/tes.bed"

  if [[ $? -ne 0 ]]; then
    echo "GTF | ERROR | Failed to generate gene catalogs"
    exit 1
  fi

  # Count entries
  GENE_COUNT=$(tail -n +2 "${WORK_DIR}/genes.tsv" | wc -l | tr -d ' ')
  TSS_COUNT=$(wc -l < "${WORK_DIR}/tss.bed" | tr -d ' ')
  TES_COUNT=$(wc -l < "${WORK_DIR}/tes.bed" | tr -d ' ')
  
  echo "GTF | PROCESS | Generated ${GENE_COUNT} gene entries"
  echo "GTF | PROCESS | Generated ${TSS_COUNT} TSS sites"
  echo "GTF | PROCESS | Generated ${TES_COUNT} TES sites"

  # Sort BED files (deterministic, coordinate-sorted)
  echo "GTF | PROCESS | Sorting BED files by coordinates..."
  LC_ALL=C sort -k1,1 -k2,2n -o "${WORK_DIR}/tss.bed" "${WORK_DIR}/tss.bed"
  LC_ALL=C sort -k1,1 -k2,2n -o "${WORK_DIR}/tes.bed" "${WORK_DIR}/tes.bed"
  echo "GTF | PROCESS | BED files sorted"

  # Move to cache
  mv -f "${WORK_DIR}/genes.tsv" "${CACHE_GENES}"
  mv -f "${WORK_DIR}/tss.bed"   "${CACHE_TSS}"
  mv -f "${WORK_DIR}/tes.bed"   "${CACHE_TES}"
  echo "GTF | PROCESS | Catalogs saved to cache"

  # Copy to output
  cp -f "${CACHE_GTF}"   "${OUT_GTF}"
  cp -f "${CACHE_GENES}" "${OUT_GENES}"
  cp -f "${CACHE_TSS}"   "${OUT_TSS}"
  cp -f "${CACHE_TES}"   "${OUT_TES}"

  ###########################################################################
  # 5) VALIDATION
  ###########################################################################

  echo "GTF | VALIDATE | Checking all output files..."
  
  ALL_VALID=1
  for FILE in "${OUT_GTF}" "${OUT_GENES}" "${OUT_TSS}" "${OUT_TES}"; do
    if [[ ! -s "${FILE}" ]]; then
      echo "GTF | ERROR | Missing or empty output: ${FILE}"
      ALL_VALID=0
    fi
  done

  if [[ ${ALL_VALID} -eq 0 ]]; then
    echo "GTF | ERROR | Not all expected outputs were created"
    exit 1
  fi

  # Report final file sizes
  echo "────────────────────────────────────────────────────────────────────────"
  echo "GTF | OUTPUT | ${OUT_GTF}: $(stat -f%z "${OUT_GTF}" 2>/dev/null || stat -c%s "${OUT_GTF}" 2>/dev/null) bytes"
  echo "GTF | OUTPUT | ${OUT_GENES}: $(stat -f%z "${OUT_GENES}" 2>/dev/null || stat -c%s "${OUT_GENES}" 2>/dev/null) bytes"
  echo "GTF | OUTPUT | ${OUT_TSS}: $(stat -f%z "${OUT_TSS}" 2>/dev/null || stat -c%s "${OUT_TSS}" 2>/dev/null) bytes"
  echo "GTF | OUTPUT | ${OUT_TES}: $(stat -f%z "${OUT_TES}" 2>/dev/null || stat -c%s "${OUT_TES}" 2>/dev/null) bytes"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "GTF | COMPLETE | genome=${ASM} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}