// ============================================================================
// fetch_and_build_index.nf — Genome Reference and Bowtie2 Index Management
// ============================================================================
//
// Purpose:
//   Fetches genome FASTA and builds Bowtie2 indices with intelligent caching
//
// Features:
//   • Smart caching - reuses existing indices across pipeline runs
//   • Local index detection - prefers ./genomes/ over downloads
//   • Automatic FASTA download from UCSC (primary + fallback URLs)
//   • Large genome support (--large-index for >4GB genomes)
//   • Deterministic outputs with storeDir for cross-run persistence
//   • Comprehensive validation of index completeness
//
// Index Reuse Strategy:
//   1. Check cache for complete index
//   2. Search local ./genomes/ directory
//   3. Download FASTA from UCSC if needed
//   4. Build new index
//
// Index Validation:
//   Ensures all 6 shards present: .1, .2, .3, .4, .rev.1, .rev.2
//   Supports both .bt2 (small) and .bt2l (large) index formats
//
// Inputs:
//   tuple(genome_id, source, fasta_in)
//     • genome_id : Assembly name (e.g., hg38, mm39, hs1)
//     • source    : Provenance info (e.g., "UCSC", "custom")
//     • fasta_in  : Path to FASTA, or "-" to force download
//
// Outputs:
//   ${params.output_dir}/00_references/${genome_id}/
//     ├── ${genome_id}.fa           — Genome FASTA
//     ├── ${genome_id}.fa.fai       — FASTA index
//     ├── ${genome_id}.genome.sizes — Chromosome sizes
//     ├── ${genome_id}.*.bt2*       — Bowtie2 index shards (in assets only)
//     └── README_index.txt          — Documentation
//
// Caching:
//   Heavy index files stored in: ${assets_dir}/references/${genome_id}/
//   Lightweight references linked to output directory
//
// ============================================================================

nextflow.enable.dsl = 2

process fetch_and_build_index {

  tag        { genome_id }
  label      'conda'
  cache      'deep'

  // Persistent storage for cross-run caching
  storeDir   { "${params.assets_dir ?: "${projectDir}/assets"}/references/${genome_id}" }
  
  // Optional: Publish lightweight reference files to results (disabled by default to save space)
  // Reference files remain available in assets/ for pipeline use
  // Enable with: --publish_references true
  publishDir(
    path: "${params.output_dir}/00_references/${genome_id}",
    mode: 'link',
    enabled: params.get('publish_references', false),
    saveAs: { filename ->
      def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
      // Exclude .bt2 index files from publishing (keep in assets only)
      if (name.endsWith('.bt2') || name.contains('.bt2l') || name.matches('.*\\.bt2(\\..*)?$')) {
        return null
      }
      return name
    }
  )

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(genome_id), val(source), path(fasta_in)

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(genome_id), val(source), path("${genome_id}.fa"), emit: ref_meta
    path "${genome_id}.fa.fai"
    path "${genome_id}.genome.sizes"
    path "${genome_id}.*.bt2*",    emit: index_files
    path "README_index.txt"

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "INDEX | START | genome=!{genome_id} | source=!{source} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  GENOME_ID="!{genome_id}"
  SOURCE="!{source}"
  THREADS=!{task.cpus}
  FORCE_REBUILD="!{params.get('force_rebuild', false) ? 'true' : 'false'}"

  # Directory structure
  CACHE_ROOT="!{params.genome_cache ?: '/tmp/genomes_cache'}"
  GENOMES_DIR="${CACHE_ROOT}/local_genomes"
  CACHE_DIR="${CACHE_ROOT}/${GENOME_ID}"
  mkdir -p "${CACHE_DIR}"

  # File paths
  FASTA_CACHE="${CACHE_DIR}/${GENOME_ID}.fa"
  FASTA_DIGEST="${CACHE_DIR}/.${GENOME_ID}.fa.sha256"
  INDEX_PREFIX="${CACHE_DIR}/${GENOME_ID}"

  # UCSC download URLs
  URL_PRIMARY="https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_ID}/bigZips/${GENOME_ID}.fa.gz"
  URL_FALLBACK="https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_ID}/bigZips/chromFa.tar.gz"

  echo "INDEX | CONFIG | Genome ID: ${GENOME_ID}"
  echo "INDEX | CONFIG | Source: ${SOURCE}"
  echo "INDEX | CONFIG | Threads: ${THREADS}"
  echo "INDEX | CONFIG | Force rebuild: ${FORCE_REBUILD}"
  echo "INDEX | CONFIG | Cache directory: ${CACHE_DIR}"
  echo "INDEX | CONFIG | Local genomes: ${GENOMES_DIR}"

  # Create temporary working directory
  TEMP_DIR=$(mktemp -d -p "${CACHE_DIR}" ".idx_${GENOME_ID}.XXXX")
  trap 'rm -rf "${TEMP_DIR}"' EXIT
  echo "INDEX | CONFIG | Temp directory: ${TEMP_DIR}"

  ###########################################################################
  # 2) VALIDATE TOOLS
  ###########################################################################

  echo "INDEX | VALIDATE | Checking required tools..."

  TOOLS_OK=1
  for TOOL in bowtie2-build samtools curl tar; do
    if command -v ${TOOL} >/dev/null 2>&1; then
      echo "INDEX | VALIDATE | ${TOOL}: $(which ${TOOL})"
    else
      echo "INDEX | ERROR | ${TOOL} not found in PATH"
      TOOLS_OK=0
    fi
  done

  if [[ ${TOOLS_OK} -eq 0 ]]; then
    echo "INDEX | ERROR | Missing required tools"
    exit 1
  fi

  ###########################################################################
  # 3) HELPER FUNCTIONS
  ###########################################################################

  # Check if complete Bowtie2 index family exists
  # Args: prefix, extension (bt2 or bt2l)
  # Returns: 0 if complete, 1 if incomplete
  has_complete_index() {
    local prefix="$1"
    local ext="$2"
    
    for shard in 1 2 3 4 rev.1 rev.2; do
      if [[ ! -s "${prefix}.${shard}.${ext}" ]]; then
        return 1
      fi
    done
    return 0
  }

  # Find local Bowtie2 index in genomes directory
  # Args: base_dir, genome_id
  # Returns: path to index prefix, or empty string
  find_local_index() {
    local base_dir="$1"
    local genome="$2"
    
    if [[ ! -d "${base_dir}" ]]; then
      echo ""
      return 0
    fi

    # Check standard locations
    for candidate in "${base_dir}/${genome}" "${base_dir}/${genome}/${genome}"; do
      if has_complete_index "${candidate}" "bt2" || has_complete_index "${candidate}" "bt2l"; then
        echo "${candidate}"
        return 0
      fi
    done

    # Search for index files with genome prefix
    shopt -s nullglob
    for ext in bt2 bt2l; do
      for index_file in "${base_dir}/${genome}"*.1."${ext}" "${base_dir}/${genome}".rev.1."${ext}"; do
        if [[ ! -s "${index_file}" ]]; then
          continue
        fi
        
        # Extract prefix
        local prefix="${index_file%.1.${ext}}"
        prefix="${prefix%.rev}"
        
        if has_complete_index "${prefix}" "${ext}"; then
          echo "${prefix}"
          shopt -u nullglob
          return 0
        fi
      done
    done
    shopt -u nullglob

    echo ""
  }

  # Stage local FASTA file if available (also copies .fai if present)
  stage_local_fasta() {
    local prefix="$1"
    
    if [[ -s "${FASTA_CACHE}" ]]; then
      return 0
    fi

    # Check various FASTA locations
    for fasta_path in \
      "${prefix}.fa" \
      "${prefix}.fa.gz" \
      "${GENOMES_DIR}/${GENOME_ID}/${GENOME_ID}.fa" \
      "${GENOMES_DIR}/${GENOME_ID}/${GENOME_ID}.fa.gz" \
      "${GENOMES_DIR}/${GENOME_ID}.fa" \
      "${GENOMES_DIR}/${GENOME_ID}.fa.gz"; do
      
      if [[ -s "${fasta_path}" ]]; then
        echo "INDEX | FASTA | Found local FASTA: ${fasta_path}"
        
        if [[ "${fasta_path}" == *.gz ]]; then
          gunzip -c "${fasta_path}" > "${FASTA_CACHE}"
        else
          cp -f "${fasta_path}" "${FASTA_CACHE}"
          
          # Also copy .fai if present
          if [[ -s "${fasta_path}.fai" ]]; then
            cp -f "${fasta_path}.fai" "${FASTA_CACHE}.fai"
            echo "INDEX | FASTA | Copied FASTA index: ${fasta_path}.fai"
          fi
        fi
        
        return 0
      fi
    done

    return 1
  }

  ###########################################################################
  # 4) CHECK FOR EXISTING INDEX (Fast Path)
  ###########################################################################

  echo "INDEX | CACHE | Checking for existing indices..."

  INDEX_READY=0

  if [[ "${FORCE_REBUILD}" == "true" ]]; then
    echo "INDEX | CACHE | Force rebuild enabled, skipping cache check"
  else
    # Check if index already in cache
    if has_complete_index "${INDEX_PREFIX}" "bt2" || has_complete_index "${INDEX_PREFIX}" "bt2l"; then
      echo "INDEX | CACHE | Found complete index in cache"
      INDEX_READY=1
    else
      # Search for local indices
      echo "INDEX | CACHE | No cached index, searching local genomes..."
      
      LOCAL_INDEX=""
      LOCAL_INDEX=$(find_local_index "${GENOMES_DIR}/${GENOME_ID}" "${GENOME_ID}")
      
      if [[ -z "${LOCAL_INDEX}" ]]; then
        LOCAL_INDEX=$(find_local_index "${GENOMES_DIR}" "${GENOME_ID}")
      fi

      if [[ -n "${LOCAL_INDEX}" ]]; then
        echo "INDEX | CACHE | Found local index: ${LOCAL_INDEX}"
        echo "INDEX | CACHE | Staging index files to cache..."
        
        # Copy all index shards
        shopt -s nullglob
        COPIED_COUNT=0
        for index_file in "${LOCAL_INDEX}".*.bt2*; do
          cp -f "${index_file}" "${CACHE_DIR}/"
          COPIED_COUNT=$((COPIED_COUNT + 1))
        done
        shopt -u nullglob
        
        echo "INDEX | CACHE | Copied ${COPIED_COUNT} index files"
        
        # Stage FASTA if available
        stage_local_fasta "${LOCAL_INDEX}"
        
        # Clean up any incomplete/zero-byte shards
        for shard in 1 2 3 4 rev.1 rev.2; do
          for ext in bt2 bt2l; do
            shard_file="${INDEX_PREFIX}.${shard}.${ext}"
            if [[ -e "${shard_file}" && ! -s "${shard_file}" ]]; then
              rm -f "${shard_file}"
            fi
          done
        done
        
        # Verify staged index is complete
        if has_complete_index "${INDEX_PREFIX}" "bt2" || has_complete_index "${INDEX_PREFIX}" "bt2l"; then
          echo "INDEX | CACHE | Local index staged successfully"
          INDEX_READY=1
        else
          echo "INDEX | CACHE | WARNING: Staged index incomplete, will build new index"
        fi
      else
        echo "INDEX | CACHE | No local indices found"
      fi
    fi
  fi

  ###########################################################################
  # 5) FETCH OR STAGE FASTA
  ###########################################################################

  if [[ ${INDEX_READY} -eq 0 ]]; then
    echo "INDEX | FASTA | Preparing genome FASTA..."

    if [[ -s "${FASTA_CACHE}" && "${FORCE_REBUILD}" != "true" ]]; then
      FASTA_SIZE=$(stat -c%s "${FASTA_CACHE}" 2>/dev/null || stat -f%z "${FASTA_CACHE}" 2>/dev/null || echo "unknown")
      echo "INDEX | FASTA | Using cached FASTA (${FASTA_SIZE} bytes)"
    else
      # Check if FASTA provided as input
      if [[ -e "!{fasta_in}" && -s "!{fasta_in}" && "!{fasta_in}" != "-" ]]; then
        echo "INDEX | FASTA | Using provided FASTA input"
        
        if [[ "!{fasta_in}" == *.gz ]]; then
          echo "INDEX | FASTA | Decompressing gzipped FASTA..."
          gunzip -c "!{fasta_in}" > "${TEMP_DIR}/${GENOME_ID}.fa"
        else
          cp -f "!{fasta_in}" "${TEMP_DIR}/${GENOME_ID}.fa"
        fi
      else
        # Download from UCSC
        echo "INDEX | FASTA | Downloading from UCSC..."
        echo "INDEX | FASTA | Primary URL: ${URL_PRIMARY}"
        
        set +e
        curl -fLsS --retry 5 --retry-delay 3 "${URL_PRIMARY}" | \
          gunzip -c > "${TEMP_DIR}/${GENOME_ID}.fa"
        DOWNLOAD_RC=$?
        set -e

        if [[ ${DOWNLOAD_RC} -ne 0 || ! -s "${TEMP_DIR}/${GENOME_ID}.fa" ]]; then
          echo "INDEX | FASTA | Primary download failed, trying fallback..."
          echo "INDEX | FASTA | Fallback URL: ${URL_FALLBACK}"
          
          curl -fLsS --retry 5 --retry-delay 3 -o "${TEMP_DIR}/chromFa.tar.gz" "${URL_FALLBACK}"
          
          mkdir -p "${TEMP_DIR}/chroms"
          tar -xzf "${TEMP_DIR}/chromFa.tar.gz" -C "${TEMP_DIR}/chroms"
          
          echo "INDEX | FASTA | Concatenating chromosome FASTAs..."
          cat ${TEMP_DIR}/chroms/*.fa > "${TEMP_DIR}/${GENOME_ID}.fa"
          
          echo "INDEX | FASTA | Fallback download successful"
        else
          echo "INDEX | FASTA | Primary download successful"
        fi
      fi

      # Validate FASTA
      if [[ ! -s "${TEMP_DIR}/${GENOME_ID}.fa" ]]; then
        echo "INDEX | ERROR | FASTA file is empty"
        exit 1
      fi

      if ! grep -q '^>' "${TEMP_DIR}/${GENOME_ID}.fa"; then
        echo "INDEX | ERROR | File does not appear to be valid FASTA format"
        exit 1
      fi

      # Count sequences
      SEQ_COUNT=$(grep -c '^>' "${TEMP_DIR}/${GENOME_ID}.fa" || echo 0)
      echo "INDEX | FASTA | Validation passed: ${SEQ_COUNT} sequences"

      # Get file size
      FASTA_SIZE=$(stat -c%s "${TEMP_DIR}/${GENOME_ID}.fa" 2>/dev/null || \
                   stat -f%z "${TEMP_DIR}/${GENOME_ID}.fa" 2>/dev/null || echo "unknown")
      echo "INDEX | FASTA | FASTA size: ${FASTA_SIZE} bytes"

      # Move to cache and create digest
      mv -f "${TEMP_DIR}/${GENOME_ID}.fa" "${FASTA_CACHE}"
      sync
      
      if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${FASTA_CACHE}" | awk '{print $1}' > "${FASTA_DIGEST}"
        echo "INDEX | FASTA | SHA256 digest created"
      fi
    fi
  fi

  ###########################################################################
  # 6) CREATE FASTA INDEX AND GENOME SIZES
  ###########################################################################

  if [[ -s "${FASTA_CACHE}" ]]; then
    echo "INDEX | FASTA | Creating FASTA index and genome sizes..."
    
    # Create .fai if missing
    if [[ ! -s "${FASTA_CACHE}.fai" ]]; then
      samtools faidx "${FASTA_CACHE}"
      echo "INDEX | FASTA | FASTA index created: ${FASTA_CACHE}.fai"
    else
      echo "INDEX | FASTA | FASTA index already exists"
    fi

    # Create genome.sizes
    cut -f1,2 "${FASTA_CACHE}.fai" > "${CACHE_DIR}/${GENOME_ID}.genome.sizes"
    
    CHR_COUNT=$(wc -l < "${CACHE_DIR}/${GENOME_ID}.genome.sizes" | tr -d ' ')
    echo "INDEX | FASTA | Genome sizes created: ${CHR_COUNT} chromosomes"
  fi

  ###########################################################################
  # 7) BUILD BOWTIE2 INDEX
  ###########################################################################

  if [[ ${INDEX_READY} -eq 0 ]]; then
    echo "INDEX | BUILD | Building Bowtie2 index..."

    # Determine if large index is needed (>4GB)
    FASTA_SIZE_BYTES=$(stat -c%s "${FASTA_CACHE}" 2>/dev/null || \
                       stat -f%z "${FASTA_CACHE}" 2>/dev/null || echo 0)
    
    LARGE_FLAG=""
    if [[ ${FASTA_SIZE_BYTES} -ge $((4 * 1024 * 1024 * 1024)) ]]; then
      LARGE_FLAG="--large-index"
      echo "INDEX | BUILD | Large genome detected (>4GB), using --large-index"
    fi

    # Check if index already complete (shouldn't happen but defensive)
    if has_complete_index "${INDEX_PREFIX}" "bt2" || has_complete_index "${INDEX_PREFIX}" "bt2l"; then
      echo "INDEX | BUILD | Index already complete, skipping build"
    else
      echo "INDEX | BUILD | Running bowtie2-build..."
      echo "INDEX | BUILD | Threads: ${THREADS}"
      echo "INDEX | BUILD | Flags: ${LARGE_FLAG:-none}"
      echo "INDEX | BUILD | This may take 30-90 minutes for large genomes..."

      # Try normal build first
      set +e
      bowtie2-build --threads ${THREADS} ${LARGE_FLAG} "${FASTA_CACHE}" "${INDEX_PREFIX}"
      BUILD_RC=$?
      set -e

      # Retry with low-memory settings if failed
      if [[ ${BUILD_RC} -ne 0 ]]; then
        echo "INDEX | BUILD | WARNING: Build failed, retrying with low-memory settings..."
        echo "INDEX | BUILD | Using single thread with conservative memory settings..."
        
        bowtie2-build --threads 1 ${LARGE_FLAG} \
                      --bmaxdivn 8 --dcv 1024 \
                      "${FASTA_CACHE}" "${INDEX_PREFIX}"
      else
        echo "INDEX | BUILD | Index build successful"
      fi

      # Verify index is complete
      if ! has_complete_index "${INDEX_PREFIX}" "bt2" && \
         ! has_complete_index "${INDEX_PREFIX}" "bt2l"; then
        echo "INDEX | ERROR | Index build incomplete, missing shards"
        ls -lh "${INDEX_PREFIX}".*.bt2* 2>/dev/null || true
        exit 1
      fi

      echo "INDEX | BUILD | Index verification passed"
    fi
  fi

  ###########################################################################
  # 8) STAGE OUTPUT FILES
  ###########################################################################

  echo "INDEX | OUTPUT | Staging output files..."

  # Remove any existing symlinks to avoid cp errors
  rm -f "${GENOME_ID}.fa" "${GENOME_ID}.fa.fai" "${GENOME_ID}.genome.sizes" 2>/dev/null || true

  # Copy reference files
  cp -f "${FASTA_CACHE}" "${GENOME_ID}.fa"
  cp -f "${FASTA_CACHE}.fai" "${GENOME_ID}.fa.fai"
  cp -f "${CACHE_DIR}/${GENOME_ID}.genome.sizes" "${GENOME_ID}.genome.sizes"

  echo "INDEX | OUTPUT | Reference files staged"

  # Copy index shards
  shopt -s nullglob
  INDEX_COUNT=0
  for index_file in "${INDEX_PREFIX}".*.bt2*; do
    cp -f "${index_file}" .
    INDEX_COUNT=$((INDEX_COUNT + 1))
  done
  shopt -u nullglob

  echo "INDEX | OUTPUT | Index shards staged: ${INDEX_COUNT} files"

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "INDEX | README | Creating documentation..."

  # Determine index type
  INDEX_TYPE="unknown"
  if has_complete_index "${INDEX_PREFIX}" "bt2"; then
    INDEX_TYPE="standard (.bt2)"
  elif has_complete_index "${INDEX_PREFIX}" "bt2l"; then
    INDEX_TYPE="large (.bt2l)"
  fi

  # Get final file sizes
  FA_SIZE=$(stat -c%s "${GENOME_ID}.fa" 2>/dev/null || \
            stat -f%z "${GENOME_ID}.fa" 2>/dev/null || echo "unknown")
  FAI_SIZE=$(stat -c%s "${GENOME_ID}.fa.fai" 2>/dev/null || \
             stat -f%z "${GENOME_ID}.fa.fai" 2>/dev/null || echo "unknown")
  SIZES_SIZE=$(stat -c%s "${GENOME_ID}.genome.sizes" 2>/dev/null || \
               stat -f%z "${GENOME_ID}.genome.sizes" 2>/dev/null || echo "unknown")

  cat > README_index.txt <<DOCEOF
================================================================================
GENOME REFERENCE AND INDEX — ${GENOME_ID}
================================================================================

REFERENCE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Genome ID:        ${GENOME_ID}
  Source:           ${SOURCE}
  Generated:        $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  
FILES
────────────────────────────────────────────────────────────────────────────
  ${GENOME_ID}.fa              — Genome FASTA (${FA_SIZE} bytes)
  ${GENOME_ID}.fa.fai          — FASTA index (${FAI_SIZE} bytes)
  ${GENOME_ID}.genome.sizes    — Chromosome sizes (${SIZES_SIZE} bytes)
  ${GENOME_ID}.*.bt2*          — Bowtie2 index shards (${INDEX_COUNT} files)

INDEX DETAILS
────────────────────────────────────────────────────────────────────────────
  Index Type:       ${INDEX_TYPE}
  Build Threads:    ${THREADS}
  Large Index:      $([ -n "${LARGE_FLAG}" ] && echo "Yes (>4GB genome)" || echo "No")
  Force Rebuild:    ${FORCE_REBUILD}

CACHING
────────────────────────────────────────────────────────────────────────────
  Cache Directory:  ${CACHE_DIR}
  Local Genomes:    ${GENOMES_DIR}
  
  The pipeline automatically caches and reuses indices across runs.
  To force a rebuild, use: --force_rebuild true

LOCAL INDEX REUSE
────────────────────────────────────────────────────────────────────────────
  The pipeline searches for existing indices in:
    1. ${CACHE_DIR}
    2. ${GENOMES_DIR}/${GENOME_ID}/
    3. ${GENOMES_DIR}/
  
  Place pre-built indices in any of these locations for automatic detection.

BOWTIE2 INDEX STRUCTURE
────────────────────────────────────────────────────────────────────────────
  Standard index (.bt2):  For genomes <4GB
    • ${GENOME_ID}.1.bt2, .2.bt2, .3.bt2, .4.bt2
    • ${GENOME_ID}.rev.1.bt2, .rev.2.bt2
  
  Large index (.bt2l):    For genomes ≥4GB
    • ${GENOME_ID}.1.bt2l, .2.bt2l, .3.bt2l, .4.bt2l
    • ${GENOME_ID}.rev.1.bt2l, .rev.2.bt2l

NOTES
────────────────────────────────────────────────────────────────────────────
  • Index files are stored in assets directory for cross-run persistence
  • Only lightweight reference files (FASTA, .fai, .sizes) are linked to results
  • Heavy index shards remain in assets to save disk space
  • SHA256 digest maintained in cache for integrity checking

================================================================================
DOCEOF

  echo "INDEX | README | Documentation created"

  ###########################################################################
  # 10) VALIDATION AND SUMMARY
  ###########################################################################

  echo "INDEX | VALIDATE | Final validation..."

  ALL_VALID=1

  # Check reference files
  for file in "${GENOME_ID}.fa" "${GENOME_ID}.fa.fai" "${GENOME_ID}.genome.sizes"; do
    if [[ ! -s "${file}" ]]; then
      echo "INDEX | ERROR | Missing or empty: ${file}"
      ALL_VALID=0
    fi
  done

  # Check index shards
  INDEX_VALID=0
  if has_complete_index "${GENOME_ID}" "bt2"; then
    INDEX_VALID=1
    INDEX_FORMAT="bt2"
  elif has_complete_index "${GENOME_ID}" "bt2l"; then
    INDEX_VALID=1
    INDEX_FORMAT="bt2l"
  fi

  if [[ ${INDEX_VALID} -eq 0 ]]; then
    echo "INDEX | ERROR | Incomplete index in output directory"
    ALL_VALID=0
  fi

  if [[ ${ALL_VALID} -eq 0 ]]; then
    echo "INDEX | ERROR | Validation failed"
    exit 1
  fi

  echo "INDEX | VALIDATE | All files validated successfully"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "INDEX | SUMMARY | Genome: ${GENOME_ID}"
  echo "INDEX | SUMMARY | FASTA: ${FA_SIZE} bytes"
  echo "INDEX | SUMMARY | Chromosomes: ${CHR_COUNT}"
  echo "INDEX | SUMMARY | Index format: ${INDEX_FORMAT}"
  echo "INDEX | SUMMARY | Index shards: ${INDEX_COUNT} files"
  echo "INDEX | SUMMARY | Cache: ${CACHE_DIR}"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "INDEX | COMPLETE | genome=${GENOME_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}