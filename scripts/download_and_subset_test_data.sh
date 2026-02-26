#!/usr/bin/env bash
# =============================================================================
# download_and_subset_test_data.sh — Download SRA FASTQs and subset to 10%
# =============================================================================
# Creates small test datasets for test_SE and test_PE. Deletes full files after
# subsetting to save space. Ships with the pipeline for easy testing.
#
# Usage:
#   ./scripts/download_and_subset_test_data.sh [SE|PE|all]
#   ./scripts/download_and_subset_test_data.sh --docker [SE|PE|all]
#
# Options:
#   --docker   Run inside TrackTx Docker image (same as pipeline)
#
# Requires: curl, gzip. Optional: fasterq-dump for NCBI fallback.
# =============================================================================

set +x 2>/dev/null || true
set -euo pipefail

# Re-exec inside Docker if --docker requested
if [[ "${1:-}" == "--docker" ]]; then
  shift
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
  DOCKER_IMAGE="ghcr.io/serhataktay/tracktx:3.0"
  echo "Running inside Docker (${DOCKER_IMAGE})..."
  exec docker run --rm --entrypoint "" -v "${PROJECT_DIR}:/data" -w /data "${DOCKER_IMAGE}" \
    bash --norc --noprofile "/data/scripts/download_and_subset_test_data.sh" "$@"
fi

FRACTION=10   # percent of reads to keep
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

# SE: GSE127844 K562
SE_SRR="SRR8669162"
SE_SAMPLE="K562_control_rep1"
SE_OUT_DIR="test_SE/test_data"

# PE: GSE154746 K562
PE_SRR="SRR12267707"
PE_SAMPLE="K562_uC0_rep1"
PE_OUT_DIR="test_PE/test_data"

# -----------------------------------------------------------------------------
# Download from ENA (direct FASTQ URLs, no sra-tools needed)
# -----------------------------------------------------------------------------
download_ena() {
  local srr="$1"
  local out_dir="$2"
  mkdir -p "${out_dir}"
  cd "${out_dir}"
  local report
  report=$(mktemp)
  if curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=fastq_ftp" -o "${report}" 2>/dev/null; then
    local ftp
    ftp=$(tail -n 1 "${report}" | cut -f2)
    rm -f "${report}"
    if [[ -z "${ftp}" || "${ftp}" == "fastq_ftp" ]]; then
      cd "${PROJECT_DIR}"
      return 1
    fi
    for url in $(echo "${ftp}" | tr ';' ' '); do
      [[ -z "${url}" ]] && continue
      local fname="${url##*/}"
      echo "  Downloading ${fname}..."
      if curl -sfL "https://${url}" -o "${fname}" 2>/dev/null; then
        echo "  Downloaded ${fname}"
      else
        curl -sfL "ftp://${url}" -o "${fname}"
      fi
    done
  else
    rm -f "${report}"
    cd "${PROJECT_DIR}"
    return 1
  fi
  cd "${PROJECT_DIR}"
  return 0
}

# -----------------------------------------------------------------------------
# Download using fasterq-dump (requires sra-tools)
# -----------------------------------------------------------------------------
download_fasterq() {
  local srr="$1"
  local out_dir="$2"
  mkdir -p "${out_dir}"
  cd "${out_dir}"
  if command -v fasterq-dump >/dev/null 2>&1; then
    echo "  Running fasterq-dump ${srr}..."
    fasterq-dump --split-files -e 4 -O . "${srr}" 2>/dev/null || true
    # Normalize names
    [[ -f "${srr}_1.fastq" ]] && mv -f "${srr}_1.fastq" "${srr}_R1.fastq"
    [[ -f "${srr}_2.fastq" ]] && mv -f "${srr}_2.fastq" "${srr}_R2.fastq"
    [[ -f "${srr}.fastq" ]] && mv -f "${srr}.fastq" "${srr}_R1.fastq"
    cd "${PROJECT_DIR}"
    return 0
  fi
  cd "${PROJECT_DIR}"
  return 1
}

# -----------------------------------------------------------------------------
# Subset FASTQ to first N% of reads (keeps complete reads only)
# -----------------------------------------------------------------------------
subset_fastq() {
  local infile="$1"
  local outfile="$2"
  local pct="${3:-10}"
  local lines total_lines target_lines
  echo "  Counting reads in ${infile}..."
  if [[ "${infile}" == *.gz ]]; then
    lines=$(gzip -dc "${infile}" | wc -l)
  else
    lines=$(wc -l < "${infile}")
  fi
  total_lines="${lines}"
  # 10% of reads = 10% of lines (4 lines per read)
  target_lines=$(( total_lines * pct / 100 ))
  # Round down to multiple of 4 for complete reads
  target_lines=$(( (target_lines / 4) * 4 ))
  [[ ${target_lines} -lt 4 ]] && target_lines=4
  echo "  Subsetting to ${target_lines} lines (~$(( target_lines / 4 )) reads)..."
  if [[ "${infile}" == *.gz ]]; then
    # head closes pipe early; gzip -dc gets SIGPIPE (141) - expected, not a real error
    gzip -dc "${infile}" | head -n "${target_lines}" | gzip -c > "${outfile}" || [[ $? -eq 141 ]]
  else
    head -n "${target_lines}" "${infile}" | gzip -c > "${outfile}" || [[ $? -eq 141 ]]
  fi
  echo "  Wrote ${outfile}"
}

# -----------------------------------------------------------------------------
# Process SE
# -----------------------------------------------------------------------------
do_se() {
  echo ""
  echo "=== test_SE: ${SE_SRR} (${SE_SAMPLE}) ==="
  mkdir -p "${SE_OUT_DIR}"
  local r1_full r1_sub
  r1_full="${SE_OUT_DIR}/${SE_SRR}_R1.fastq.gz"
  r1_sub="${SE_OUT_DIR}/${SE_SAMPLE}_R1.fastq.gz"
  if [[ -f "${r1_sub}" ]]; then
    echo "  Subset already exists: ${r1_sub}"
    return 0
  fi
  if [[ ! -f "${r1_full}" ]]; then
    # Try ENA first (often has .gz directly)
    if ! download_ena "${SE_SRR}" "${SE_OUT_DIR}"; then
      echo "  ENA failed, trying fasterq-dump..."
      if ! download_fasterq "${SE_SRR}" "${SE_OUT_DIR}"; then
        echo "ERROR: Download failed. Install sra-tools (conda install -c bioconda sra-tools)"
        return 1
      fi
    fi
    # Find R1 (ENA may use .fastq.gz, _1.fastq.gz; fasterq-dump uses _R1.fastq)
    if [[ -f "${SE_OUT_DIR}/${SE_SRR}.fastq.gz" ]]; then
      r1_full="${SE_OUT_DIR}/${SE_SRR}.fastq.gz"
    elif [[ -f "${SE_OUT_DIR}/${SE_SRR}_1.fastq.gz" ]]; then
      r1_full="${SE_OUT_DIR}/${SE_SRR}_1.fastq.gz"
    elif [[ -f "${SE_OUT_DIR}/${SE_SRR}_R1.fastq" ]]; then
      gzip -c "${SE_OUT_DIR}/${SE_SRR}_R1.fastq" > "${r1_full}"
      rm -f "${SE_OUT_DIR}/${SE_SRR}_R1.fastq"
    fi
  fi
  subset_fastq "${r1_full}" "${r1_sub}" "${FRACTION}"
  rm -f "${r1_full}"
  echo "  Done: ${r1_sub} (removed full file to save space)"
}

# -----------------------------------------------------------------------------
# Process PE
# -----------------------------------------------------------------------------
do_pe() {
  echo ""
  echo "=== test_PE: ${PE_SRR} (${PE_SAMPLE}) ==="
  mkdir -p "${PE_OUT_DIR}"
  local r1_full r2_full r1_sub r2_sub
  r1_full="${PE_OUT_DIR}/${PE_SRR}_R1.fastq.gz"
  r2_full="${PE_OUT_DIR}/${PE_SRR}_R2.fastq.gz"
  # ENA uses _1/_2, fasterq-dump uses _R1/_R2; resolve to actual files
  [[ -f "${PE_OUT_DIR}/${PE_SRR}_1.fastq.gz" ]] && r1_full="${PE_OUT_DIR}/${PE_SRR}_1.fastq.gz"
  [[ -f "${PE_OUT_DIR}/${PE_SRR}_2.fastq.gz" ]] && r2_full="${PE_OUT_DIR}/${PE_SRR}_2.fastq.gz"
  r1_sub="${PE_OUT_DIR}/${PE_SAMPLE}_R1.fastq.gz"
  r2_sub="${PE_OUT_DIR}/${PE_SAMPLE}_R2.fastq.gz"
  if [[ -f "${r1_sub}" && -f "${r2_sub}" ]]; then
    echo "  Subset already exists: ${r1_sub}, ${r2_sub}"
    return 0
  fi
  if [[ ! -f "${r1_full}" || ! -f "${r2_full}" ]]; then
    if ! download_ena "${PE_SRR}" "${PE_OUT_DIR}"; then
      echo "  ENA failed, trying fasterq-dump..."
      if ! download_fasterq "${PE_SRR}" "${PE_OUT_DIR}"; then
        echo "ERROR: Download failed. Install sra-tools (conda install -c bioconda sra-tools)"
        return 1
      fi
    fi
    # Resolve R1/R2 (ENA: _1.fastq.gz, _2.fastq.gz; fasterq-dump: _R1.fastq, _R2.fastq)
    if [[ -f "${PE_OUT_DIR}/${PE_SRR}_1.fastq.gz" ]]; then
      r1_full="${PE_OUT_DIR}/${PE_SRR}_1.fastq.gz"
    elif [[ -f "${PE_OUT_DIR}/${PE_SRR}_R1.fastq" ]]; then
      gzip -c "${PE_OUT_DIR}/${PE_SRR}_R1.fastq" > "${r1_full}"
      rm -f "${PE_OUT_DIR}/${PE_SRR}_R1.fastq"
    fi
    if [[ -f "${PE_OUT_DIR}/${PE_SRR}_2.fastq.gz" ]]; then
      r2_full="${PE_OUT_DIR}/${PE_SRR}_2.fastq.gz"
    elif [[ -f "${PE_OUT_DIR}/${PE_SRR}_R2.fastq" ]]; then
      gzip -c "${PE_OUT_DIR}/${PE_SRR}_R2.fastq" > "${r2_full}"
      rm -f "${PE_OUT_DIR}/${PE_SRR}_R2.fastq"
    fi
  fi
  if [[ ! -f "${r1_full}" ]]; then
    echo "ERROR: R1 not found: ${r1_full}"
    return 1
  fi
  if [[ ! -f "${r2_full}" ]]; then
    echo "ERROR: R2 not found: ${r2_full} (PE requires both reads)"
    return 1
  fi
  subset_fastq "${r1_full}" "${r1_sub}" "${FRACTION}"
  subset_fastq "${r2_full}" "${r2_sub}" "${FRACTION}"
  rm -f "${r1_full}" "${r2_full}"
  echo "  Done: ${r1_sub}, ${r2_sub} (removed full files to save space)"
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
MODE="${1:-all}"
echo "TrackTx test data: download + subset to ${FRACTION}%"
case "${MODE}" in
  SE|se)  do_se ;;
  PE|pe)  do_pe ;;
  all)    do_se; do_pe ;;
  -h|--help)
    echo "Usage: $0 [--docker] [SE|PE|all]"
    echo "  --docker  Run inside TrackTx Docker image"
    echo "  SE        Single-end test data only"
    echo "  PE        Paired-end test data only"
    echo "  all       Both (default)"
    exit 0 ;;
  *)
    echo "Usage: $0 [--docker] [SE|PE|all]"
    exit 1 ;;
esac
echo ""
echo "Done. Run pipeline with:"
echo "  SE: ./run_pipeline.sh --samplesheet test_SE/samplesheet_SE.csv --params-file test_SE/params_SE.yaml --output_dir ./results_test_SE"
echo "  PE: ./run_pipeline.sh --samplesheet test_PE/samplesheet_PE.csv --params-file test_PE/params_PE.yaml --output_dir ./results_test_PE"
