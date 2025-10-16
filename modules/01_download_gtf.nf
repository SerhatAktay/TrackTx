// ============================================================================
// download_gtf.nf — GTF fetch + gene catalog + functional regions
// ----------------------------------------------------------------------------
// Overview
//   • Uses params.* to obtain a reference GTF (custom path/URL or UCSC)
//   • Caches artifacts (GTF + derived tables) under params.genome_cache
//   • Creates gene catalog and functional region coordinate files
//   • Persists deterministic outputs via storeDir; publishes links under results/
//   • Derived BED6 files are sorted (chr, start asc)
//
// Inputs
//   (none — driven by params)
//   Required params:
//     params.reference_genome : hs1|hg38|mm10|mm39|dm6|TAIR10|other
//     params.genome_cache     : cache dir for fetched/derived files
//     params.assets_dir       : canonical artifact store (for storeDir)
//   When reference_genome == "other":
//     params.gtf_path | params.gtf_url
//
// Outputs (publishDir):
//   ${params.output_dir}/00_references/${params.reference_genome}/
//     ├── ${reference_genome}.gtf
//     ├── ${reference_genome}.genes.tsv  (gene_id, gene_name, chr, strand, start, end, tss, tes, biotype)
//     ├── ${reference_genome}.tss.bed    (BED6; 1bp TSS per gene)
//     ├── ${reference_genome}.tes.bed    (BED6; 1bp TES per gene)
//     └── ${reference_genome}.functional_regions.tsv  (all functional region coordinates)
// ============================================================================

nextflow.enable.dsl = 2

process download_gtf {

  // ── Meta / resources ─────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  tag        { params.reference_genome }
  label      'conda'
  cache      'lenient'

  // Use the main tracktx conda environment
  conda 'envs/tracktx.yaml'

  // Persist artifacts; publish lightweight links
  storeDir   "${params.assets_dir ?: "${projectDir}/assets"}/annotation/${params.reference_genome}"
  publishDir "${params.output_dir}/00_references/${params.reference_genome}", mode: 'link', overwrite: true

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    // none

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    path("${params.reference_genome}.gtf")
    path("${params.reference_genome}.genes.tsv"), emit: genes
    path("${params.reference_genome}.tss.bed"),  emit: tss
    path("${params.reference_genome}.tes.bed"),  emit: tes
    path("${params.reference_genome}.functional_regions.tsv"), emit: functional_regions

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  trap 'echo "ERROR [gtf] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  export LC_ALL=C
  
  # Debug: show which Python is being used
  echo "DEBUG [gtf] Using Python: $(which python3)"
  echo "DEBUG [gtf] Python version: $(python3 --version)"

  # ── Bindings / locals ────────────────────────────────────────────────────
  ASM='!{params.reference_genome}'
  CACHE_DIR='!{ params.genome_cache ?: "/tmp/genomes_cache" }'
  mkdir -p "${CACHE_DIR}"

  OUT_GTF="${ASM}.gtf"
  OUT_GENES="${ASM}.genes.tsv"
  OUT_TSS="${ASM}.tss.bed"
  OUT_TES="${ASM}.tes.bed"
  OUT_FUNC_REGIONS="${ASM}.functional_regions.tsv"

  CACHE_GTF="${CACHE_DIR}/${ASM}.gtf"
  CACHE_GENES="${CACHE_DIR}/${ASM}.genes.tsv"
  CACHE_TSS="${CACHE_DIR}/${ASM}.tss.bed"
  CACHE_TES="${CACHE_DIR}/${ASM}.tes.bed"
  CACHE_FUNC_REGIONS="${CACHE_DIR}/${ASM}.functional_regions.tsv"

  CUSTOM_PATH='!{params.gtf_path ?: ""}'
  CUSTOM_URL='!{params.gtf_url  ?: ""}'

  RSYNC_URL="rsync://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/bigZips/genes/${ASM}.ncbiRefSeq.gtf.gz"
  HTTPS_URL="https://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/bigZips/genes/${ASM}.ncbiRefSeq.gtf.gz"

  # Local GTF preference (mirrors index module's local reuse)
  # If a GTF exists under ./genomes/<ASM>.gtf or ./genomes/<ASM>/<ASM>.gtf, use it.
  stage_local_gtf_if_any(){
    local asm="$1"
    # Look for local genomes in standard locations
    local g1="${CACHE_DIR}/local_genomes/${asm}.gtf"
    local g2="${CACHE_DIR}/local_genomes/${asm}/${asm}.gtf"
    if [[ -s "$g1" ]]; then
      cp -f "$g1" "${CACHE_GTF}"; return 0
    elif [[ -s "$g2" ]]; then
      cp -f "$g2" "${CACHE_GTF}"; return 0
    fi
    return 1
  }

  echo "INFO  [gtf] ▶ assembly=${ASM} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  command -v python3 >/dev/null 2>&1 || { echo "ERROR python3 not found" >&2; exit 1; }

  # ─────────────────────────────────────────────────────────────────────────
  # 1) Fast path: reuse complete cache
  # ─────────────────────────────────────────────────────────────────────────
  if [[ -s "${CACHE_GTF}" && -s "${CACHE_GENES}" && -s "${CACHE_TSS}" && -s "${CACHE_TES}" && -s "${CACHE_FUNC_REGIONS}" ]]; then
    echo "INFO  [gtf] using cached GTF + derived catalogs + functional regions"
    cp -f "${CACHE_GTF}"   "${OUT_GTF}"
    cp -f "${CACHE_GENES}" "${OUT_GENES}"
    cp -f "${CACHE_TSS}"   "${OUT_TSS}"
    cp -f "${CACHE_TES}"   "${OUT_TES}"
    cp -f "${CACHE_FUNC_REGIONS}" "${OUT_FUNC_REGIONS}"
    exit 0
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 2) Ensure GTF present in cache (custom → rsync → https)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ -s "${CACHE_GTF}" ]]; then
    echo "INFO  [gtf] cached GTF found; rebuilding derived tables"
    cp -f "${CACHE_GTF}" "${OUT_GTF}"
  else
    tmpd="$(mktemp -d "${CACHE_DIR}/.gtf_${ASM}.XXXXXX")"; trap 'rm -rf "${tmpd}"' EXIT
    GTF_TMP="${tmpd}/${ASM}.gtf"
    fetch_ok=0

    if [[ "${ASM}" == "other" && -n "${CUSTOM_PATH}" && -e "${CUSTOM_PATH}" ]]; then
      echo "INFO  [gtf] copying custom GTF from path"
      if [[ "${CUSTOM_PATH}" == *.gz ]]; then gzip -cd -- "${CUSTOM_PATH}" > "${GTF_TMP}"; else cp -f "${CUSTOM_PATH}" "${GTF_TMP}"; fi
      fetch_ok=1
    elif [[ "${ASM}" == "other" && -n "${CUSTOM_URL}" ]]; then
      echo "INFO  [gtf] fetching custom GTF from URL"
      if [[ "${CUSTOM_URL}" == *.gz ]]; then
        curl -fsSL --retry 3 --retry-delay 4 "${CUSTOM_URL}" | gzip -cd > "${GTF_TMP}" && fetch_ok=1 || true
      else
        curl -fsSL --retry 3 --retry-delay 4 "${CUSTOM_URL}" > "${GTF_TMP}" && fetch_ok=1 || true
      fi
    else
      if stage_local_gtf_if_any "${ASM}"; then
        echo "INFO  [gtf] using local genomes/${ASM}.gtf"
        cp -f "${CACHE_GTF}" "${OUT_GTF}"
        fetch_ok=1
      fi
      if [[ "${fetch_ok}" -ne 1 ]]; then
      if command -v rsync >/dev/null 2>&1; then
        echo "INFO  [gtf] UCSC rsync…"
        rsync --quiet "${RSYNC_URL}" - 2>/dev/null | gunzip -c > "${GTF_TMP}" && fetch_ok=1 || true
      fi
      if [[ "${fetch_ok}" -ne 1 ]]; then
        echo "WARN  [gtf] rsync unavailable/failed — trying HTTPS"
        curl -fsSL --retry 3 --retry-delay 4 "${HTTPS_URL}" | gunzip -c > "${GTF_TMP}" && fetch_ok=1 || true
      fi
      fi
    fi

    [[ "${fetch_ok}" -eq 1 && -s "${GTF_TMP}" ]] || { echo "ERROR [gtf] unable to obtain GTF for ${ASM}" >&2; exit 2; }

    # Minimal validation (non-fatal)
    if ! awk 'BEGIN{FS="\\t"} NF>=9{ok=1; exit} END{exit (ok?0:1)}' "${GTF_TMP}"; then
      echo "WARN  [gtf] GTF may not be strictly 9 columns; proceeding" >&2
    fi

    if [[ -s "${GTF_TMP}" ]]; then
      mv -f "${GTF_TMP}" "${CACHE_GTF}"
      cp -f "${CACHE_GTF}" "${OUT_GTF}"
    fi
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 3) Build catalogs (deterministic; BEDs sorted)
  #    - Includes ALL gene types (maintains scientific completeness)
  #    - Calculates single TSS/TES coordinates from transcript boundaries
  #    - Improved GTF parsing with better biotype handling
  # ─────────────────────────────────────────────────────────────────────────
  tmpw="$(mktemp -d "${CACHE_DIR}/.build_${ASM}.XXXXXX")"; trap 'rm -rf "${tmpw}"' EXIT

  # Use container-safe path to Python script
  BIN_SCRIPT="!{projectDir}/bin/gtf_to_catalog.py"
  if [[ ! -f "$BIN_SCRIPT" ]]; then
    BIN_SCRIPT="gtf_to_catalog.py"  # Fall back to PATH
  fi
  python3 "$BIN_SCRIPT" "${CACHE_GTF}" "${tmpw}/genes.tsv" "${tmpw}/tss.bed" "${tmpw}/tes.bed"

  LC_ALL=C sort -k1,1 -k2,2n -o "${tmpw}/tss.bed" "${tmpw}/tss.bed"
  LC_ALL=C sort -k1,1 -k2,2n -o "${tmpw}/tes.bed" "${tmpw}/tes.bed"

  mv -f "${tmpw}/genes.tsv" "${CACHE_GENES}"
  mv -f "${tmpw}/tss.bed"   "${CACHE_TSS}"
  mv -f "${tmpw}/tes.bed"   "${CACHE_TES}"

  # ─────────────────────────────────────────────────────────────────────────
  # 4) Create functional regions based on user parameters
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  [gtf] Creating functional regions based on parameters"
  
  # Use container-safe path to Python script
  FUNC_REGIONS_SCRIPT="!{projectDir}/bin/create_functional_regions.py"
  if [[ ! -f "$FUNC_REGIONS_SCRIPT" ]]; then
    FUNC_REGIONS_SCRIPT="create_functional_regions.py"  # Fall back to PATH
  fi
  
  # Get functional region parameters (with defaults matching old R script)
  PROM_UP="!{ (params.functional_regions?.prom_up ?: 250) as int }"
  PROM_DOWN="!{ (params.functional_regions?.prom_down ?: 250) as int }"
  DIV_INNER="!{ (params.functional_regions?.div_inner ?: 251) as int }"
  DIV_OUTER="!{ (params.functional_regions?.div_outer ?: 750) as int }"
  CPS_OFFSET="!{ (params.functional_regions?.cps_offset ?: 500) as int }"
  TW_LENGTH="!{ (params.functional_regions?.tw_length ?: 10000) as int }"
  
  python3 "$FUNC_REGIONS_SCRIPT" \
    --genes "${CACHE_GENES}" \
    --prom-up "$PROM_UP" \
    --prom-down "$PROM_DOWN" \
    --div-inner "$DIV_INNER" \
    --div-outer "$DIV_OUTER" \
    --cps-offset "$CPS_OFFSET" \
    --tw-length "$TW_LENGTH" \
    --output "${CACHE_FUNC_REGIONS}"

  cp -f "${CACHE_GTF}"   "${OUT_GTF}"
  cp -f "${CACHE_GENES}" "${OUT_GENES}"
  cp -f "${CACHE_TSS}"   "${OUT_TSS}"
  cp -f "${CACHE_TES}"   "${OUT_TES}"
  cp -f "${CACHE_FUNC_REGIONS}" "${OUT_FUNC_REGIONS}"

  # ─────────────────────────────────────────────────────────────────────────
  # 5) Guards
  # ─────────────────────────────────────────────────────────────────────────
  [[ -s "${OUT_GTF}" && -s "${OUT_GENES}" && -s "${OUT_TSS}" && -s "${OUT_TES}" && -s "${OUT_FUNC_REGIONS}" ]] \
    || { echo "ERROR [gtf] expected outputs missing" >&2; exit 3; }

  echo "INFO  [gtf] ✔ made: ${CACHE_GTF}, ${CACHE_GENES}, ${CACHE_TSS}, ${CACHE_TES}, ${CACHE_FUNC_REGIONS} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
