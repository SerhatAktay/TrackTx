// ============================================================================
// calculate_pol2_metrics.nf — Pol-II Metrics Calculation
// ============================================================================
//
// Purpose:
//   Calculates RNA Polymerase II metrics from alignments and coverage tracks
//
// Features:
//   • Density metrics: Signal per functional region from normalized tracks
//   • Gene-level metrics: TSS and gene body coverage from BAM
//   • Pausing Index: TSS/body ratio indicating promoter-proximal pausing
//   • Quality control: Per-sample QC JSON
//
// Two Complementary Approaches:
//   1. Density (from normalized bedGraphs):
//      - Uses siCPM (preferred) or CPM tracks
//      - Sums |pos| + |neg| signal per functional region
//      - Fast, region-level quantification
//
//   2. Gene metrics (from BAM):
//      - MAPQ-filtered alignments
//      - Optional duplicate removal
//      - TSS window and gene body coverage
//      - Pausing index calculation
//
// Pausing Index (PI):
//   PI = (TSS_density / TSS_width) / (Body_density / Body_width)
//   
//   High PI (>1.5): Strong promoter-proximal pausing
//   Low PI (<1.0):  Productive elongation
//
// Inputs:
//   tuple(sample_id, bam, functional_bed,
//         pos3_cpm_bg, neg3_cpm_bg, pos3_sicpm_bg, neg3_sicpm_bg,
//         condition, timepoint, replicate)
//   path(gtf) : Gene annotations
//
// Outputs:
//   ${params.output_dir}/08_pol2_metrics/${sample_id}/
//     ├── pol2_gene_metrics.tsv    — Per-gene TSS and body metrics
//     ├── pausing_index.tsv        — Per-gene pausing indices
//     ├── pol2_density.tsv         — Signal per functional region
//     ├── pol2_qc.json             — QC summary (optional)
//     ├── README_pol2_metrics.txt  — Documentation
//     └── pol2_metrics.log         — Processing log
//
// Parameters (params.pol2.*):
//   mapq               : Minimum MAPQ (default: 10)
//   dedup              : Remove duplicates (default: true)
//   tss_win            : TSS window ±bp (default: 50)
//   body_offset_min    : Min body start offset (default: 2000 bp)
//   body_offset_frac   : Body start as fraction (default: 0.10)
//   feature_types      : GTF features to use (default: "gene,transcript")
//   fail_if_no_genes   : Fail on empty gene list (default: false)
//
// ============================================================================

nextflow.enable.dsl = 2

process calculate_pol2_metrics {

  tag        { sid }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/08_pol2_metrics/${sid}",
             mode: 'copy',
             overwrite: true

  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sid),
          path(in_bam),
          path(func_bed),
          path(pos3_cpm_bg), path(neg3_cpm_bg),
          path(pos3_sicpm_bg), path(neg3_sicpm_bg),
          val(cond), val(tp), val(rep)
    path gtf

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sid), path('pol2_gene_metrics.tsv'), 
          val(cond), val(tp), val(rep), emit: genes
    tuple val(sid), path('pausing_index.tsv'),
          val(cond), val(tp), val(rep), emit: pausing
    tuple val(sid), path('pol2_density.tsv'),
          val(cond), val(tp), val(rep), emit: density
    path 'pol2_qc.json', optional: true, emit: qc
    path 'README_pol2_metrics.txt'
    path 'pol2_metrics.log', emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Limit threading for linear algebra libraries
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export BLIS_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1

  # Redirect all output to log file
  exec > >(tee -a pol2_metrics.log) 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "POL2 | START | sample=!{sid} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sid}"
  CONDITION="!{cond}"
  TIMEPOINT="!{tp}"
  REPLICATE="!{rep}"
  THREADS=!{task.cpus}

  # Input files
  IN_BAM="!{in_bam}"
  FUNC_BED="!{func_bed}"
  GTF_FILE="!{gtf}"
  CALC_SCRIPT="!{projectDir}/bin/calculate_pol2_metrics.py"

  # Coverage tracks
  POS_CPM="!{pos3_cpm_bg}"
  NEG_CPM="!{neg3_cpm_bg}"
  POS_SICPM="!{pos3_sicpm_bg}"
  NEG_SICPM="!{neg3_sicpm_bg}"

  # Parameters
  MAPQ=!{params.pol2?.mapq ?: 10}
  DEDUP_ENABLED=$([[ "!{params.pol2?.dedup ?: true}" == "false" ]] && echo 0 || echo 1)
  TSS_WIN=!{params.pol2?.tss_win ?: 50}
  BODY_OFFSET_MIN=!{params.pol2?.body_offset_min ?: 2000}
  BODY_OFFSET_FRAC=!{params.pol2?.body_offset_frac ?: 0.10}
  FEATURE_TYPES="!{params.pol2?.feature_types ?: 'gene,transcript'}"
  FAIL_IF_NO_GENES=$([[ "!{params.pol2?.fail_if_no_genes}" == "true" ]] && echo 1 || echo 0)

  echo "POL2 | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "POL2 | CONFIG | Condition: ${CONDITION}"
  echo "POL2 | CONFIG | Timepoint: ${TIMEPOINT}"
  echo "POL2 | CONFIG | Replicate: ${REPLICATE}"
  echo "POL2 | CONFIG | Threads: ${THREADS}"
  echo ""
  echo "POL2 | CONFIG | Input Files:"
  echo "POL2 | CONFIG |   BAM: $(basename ${IN_BAM})"
  echo "POL2 | CONFIG |   Functional regions: $(basename ${FUNC_BED})"
  echo "POL2 | CONFIG |   GTF: $(basename ${GTF_FILE})"
  echo ""
  echo "POL2 | CONFIG | Coverage Tracks:"
  echo "POL2 | CONFIG |   CPM: $(basename ${POS_CPM}), $(basename ${NEG_CPM})"
  echo "POL2 | CONFIG |   siCPM: $(basename ${POS_SICPM}), $(basename ${NEG_SICPM})"
  echo ""
  echo "POL2 | CONFIG | BAM Filtering:"
  echo "POL2 | CONFIG |   MAPQ threshold: ${MAPQ}"
  echo "POL2 | CONFIG |   Remove duplicates: $([ ${DEDUP_ENABLED} -eq 1 ] && echo "yes" || echo "no")"
  echo ""
  echo "POL2 | CONFIG | Gene Metrics:"
  echo "POL2 | CONFIG |   TSS window: ±${TSS_WIN} bp"
  echo "POL2 | CONFIG |   Body offset min: ${BODY_OFFSET_MIN} bp"
  echo "POL2 | CONFIG |   Body offset fraction: ${BODY_OFFSET_FRAC}"
  echo "POL2 | CONFIG |   Feature types: ${FEATURE_TYPES}"
  echo "POL2 | CONFIG |   Fail if no genes: $([ ${FAIL_IF_NO_GENES} -eq 1 ] && echo "yes" || echo "no")"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "POL2 | VALIDATE | Checking input files..."

  VALIDATION_OK=1

  # Check Python script
  if [[ ! -f "${CALC_SCRIPT}" ]]; then
    echo "POL2 | ERROR | Python script not found: ${CALC_SCRIPT}"
    VALIDATION_OK=0
  else
    echo "POL2 | VALIDATE | Python script: ${CALC_SCRIPT}"
  fi

  # Check BAM
  if [[ ! -s "${IN_BAM}" ]]; then
    echo "POL2 | ERROR | BAM file missing or empty: ${IN_BAM}"
    VALIDATION_OK=0
  else
    BAM_SIZE=$(stat -c%s "${IN_BAM}" 2>/dev/null || stat -f%z "${IN_BAM}" 2>/dev/null || echo "unknown")
    echo "POL2 | VALIDATE | BAM: ${BAM_SIZE} bytes"
  fi

  # Check GTF
  if [[ ! -e "${GTF_FILE}" ]]; then
    echo "POL2 | ERROR | GTF file missing: ${GTF_FILE}"
    VALIDATION_OK=0
  else
    GTF_SIZE=$(stat -c%s "${GTF_FILE}" 2>/dev/null || stat -f%z "${GTF_FILE}" 2>/dev/null || echo "unknown")
    GTF_LINES=$(wc -l < "${GTF_FILE}" 2>/dev/null | tr -d ' ' || echo 0)
    echo "POL2 | VALIDATE | GTF: ${GTF_SIZE} bytes (${GTF_LINES} lines)"
  fi

  # Check CPM tracks (required)
  for TRACK in "${POS_CPM}" "${NEG_CPM}"; do
    if [[ ! -e "${TRACK}" ]]; then
      echo "POL2 | ERROR | Required CPM track missing: ${TRACK}"
      VALIDATION_OK=0
    else
      TRACK_SIZE=$(stat -c%s "${TRACK}" 2>/dev/null || stat -f%z "${TRACK}" 2>/dev/null || echo "unknown")
      echo "POL2 | VALIDATE | $(basename ${TRACK}): ${TRACK_SIZE} bytes"
    fi
  done

  # Check siCPM tracks (optional)
  if [[ -s "${POS_SICPM}" && -s "${NEG_SICPM}" ]]; then
    SICPM_SIZE=$(stat -c%s "${POS_SICPM}" 2>/dev/null || stat -f%z "${POS_SICPM}" 2>/dev/null || echo "unknown")
    echo "POL2 | VALIDATE | siCPM tracks available: ${SICPM_SIZE} bytes"
    SICPM_AVAILABLE=1
  else
    echo "POL2 | VALIDATE | siCPM tracks not available, will use CPM"
    SICPM_AVAILABLE=0
  fi

  # Check functional regions (optional but expected)
  if [[ "${FUNC_BED}" != "-" && -s "${FUNC_BED}" ]]; then
    FUNC_SIZE=$(stat -c%s "${FUNC_BED}" 2>/dev/null || stat -f%z "${FUNC_BED}" 2>/dev/null || echo "unknown")
    FUNC_COUNT=$(grep -v '^#' "${FUNC_BED}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "POL2 | VALIDATE | Functional regions: ${FUNC_COUNT} regions (${FUNC_SIZE} bytes)"
  else
    echo "POL2 | VALIDATE | WARNING: No functional regions provided"
  fi

  # Validate tools
  for TOOL in samtools bedtools python3 awk; do
    if command -v ${TOOL} >/dev/null 2>&1; then
      echo "POL2 | VALIDATE | ${TOOL}: $(which ${TOOL})"
    else
      echo "POL2 | ERROR | Required tool not found: ${TOOL}"
      VALIDATION_OK=0
    fi
  done

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "POL2 | ERROR | Validation failed"
    exit 1
  fi

  ###########################################################################
  # 3) DENSITY METRICS FROM NORMALIZED TRACKS
  ###########################################################################

  echo "POL2 | DENSITY | Calculating density metrics from normalized tracks..."

  # Choose between siCPM (preferred) and CPM
  select_track() {
    local sicpm="$1"
    local cpm="$2"
    
    if [[ -s "${sicpm}" ]]; then
      echo "${sicpm}"
      return 0
    else
      echo "${cpm}"
      return 0
    fi
  }

  POS_TRACK=$(select_track "${POS_SICPM}" "${POS_CPM}")
  NEG_TRACK=$(select_track "${NEG_SICPM}" "${NEG_CPM}")

  if [[ "${POS_TRACK}" == "${POS_SICPM}" ]]; then
    NORM_METHOD="siCPM"
    echo "POL2 | DENSITY | Using siCPM normalization (spike-in)"
  else
    NORM_METHOD="CPM"
    echo "POL2 | DENSITY | Using CPM normalization (standard)"
  fi

  # Helper to read bedGraph (handles gzip)
  read_bedgraph() {
    local file="$1"
    
    if [[ "${file}" == *.gz ]]; then
      gzip -cd "${file}"
    else
      cat "${file}"
    fi
  }

  # Helper to clean and take absolute values
  clean_and_abs() {
    local input="$1"
    local output="$2"
    
    if [[ ! -s "${input}" ]]; then
      : > "${output}"
      return 0
    fi
    
    echo "POL2 | DENSITY | Cleaning: $(basename ${input})"
    
    read_bedgraph "${input}" | \
      awk 'BEGIN{OFS="\t"}
           /^#/ || /^track/ || /^browser/ {next}
           (NF>=4) {
             chr=$1
             start=$2+0
             end=$3+0
             val=$4+0
             if (end > start) {
               if (val < 0) val = -val
               print chr, start, end, val
             }
           }' | \
      LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${output}"
    
    LINE_COUNT=$(wc -l < "${output}" | tr -d ' ')
    echo "POL2 | DENSITY | Cleaned: ${LINE_COUNT} intervals"
  }

  # Clean and take absolute values of both strands
  clean_and_abs "${POS_TRACK}" "pos.abs.bedgraph"
  clean_and_abs "${NEG_TRACK}" "neg.abs.bedgraph"

  # Merge positive and negative strands (|pos| + |neg|)
  echo "POL2 | DENSITY | Merging strands..."
  
  cat pos.abs.bedgraph neg.abs.bedgraph | \
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n | \
    bedtools merge -i - -c 4 -o sum > combined.norm.bedgraph || \
    : > combined.norm.bedgraph

  COMBINED_LINES=$(wc -l < combined.norm.bedgraph | tr -d ' ')
  COMBINED_SIZE=$(stat -c%s combined.norm.bedgraph 2>/dev/null || stat -f%z combined.norm.bedgraph 2>/dev/null || echo "unknown")
  echo "POL2 | DENSITY | Combined track: ${COMBINED_LINES} intervals (${COMBINED_SIZE} bytes)"

  # Map signal to functional regions
  if [[ "${FUNC_BED}" != "-" && -s "${FUNC_BED}" ]]; then
    echo "POL2 | DENSITY | Mapping signal to functional regions..."
    
    # Clean functional regions BED
    awk 'BEGIN{OFS="\t"}
         !/^track/ && !/^browser/ && !/^#/ && (NF>=3) {
           print $1, $2, $3, (NF>=4 ? $4 : "."), (NF>=5 ? $5 : "0"), (NF>=6 ? $6 : ".")
         }' "${FUNC_BED}" | \
      LC_ALL=C sort -k1,1 -k2,2n -k3,3n > functional_regions.sorted.bed
    
    FUNC_SORTED=$(wc -l < functional_regions.sorted.bed | tr -d ' ')
    echo "POL2 | DENSITY | Sorted functional regions: ${FUNC_SORTED}"
    
    # Create header
    echo -e "chr\tstart\tend\tname\tsignal\tnorm_method" > pol2_density.tsv
    
    # Map signal using bedtools
    bedtools map \
      -a functional_regions.sorted.bed \
      -b combined.norm.bedgraph \
      -c 4 \
      -o sum \
      -null 0 | \
      awk -v OFS='\t' -v METHOD="${NORM_METHOD}" '{
        print $1, $2, $3, ($4 != "." ? $4 : "."), ($NF != "." ? $NF : 0), METHOD
      }' >> pol2_density.tsv
    
    DENSITY_LINES=$(tail -n +2 pol2_density.tsv | wc -l | tr -d ' ')
    echo "POL2 | DENSITY | Density table: ${DENSITY_LINES} regions"
  else
    echo "POL2 | DENSITY | No functional regions, creating header-only file"
    echo -e "chr\tstart\tend\tname\tsignal\tnorm_method" > pol2_density.tsv
  fi

  ###########################################################################
  # 4) PREPARE BAM FOR GENE METRICS
  ###########################################################################

  echo "POL2 | BAM | Preparing BAM for gene-level metrics..."

  # Check if BAM is coordinate sorted
  SO_COORD=0
  if samtools view -H "${IN_BAM}" | \
     awk '/^@HD/ && /SO:coordinate/ {ok=1} END{exit ok?0:1}'; then
    SO_COORD=1
    echo "POL2 | BAM | BAM is coordinate-sorted"
  else
    echo "POL2 | BAM | BAM is not coordinate-sorted, will sort"
  fi

  # Build filtering flags
  FILTER_FLAGS="-F 0x4"  # Exclude unmapped
  if [[ ${DEDUP_ENABLED} -eq 1 ]]; then
    FILTER_FLAGS="${FILTER_FLAGS} -F 0x400"  # Exclude duplicates
    echo "POL2 | BAM | Will exclude duplicates"
  else
    echo "POL2 | BAM | Will retain duplicates"
  fi

  # Filter BAM
  BAM_START=$(date +%s)
  
  if [[ ${SO_COORD} -eq 1 ]]; then
    echo "POL2 | BAM | Filtering BAM (MAPQ≥${MAPQ})..."
    samtools view \
      -@ ${THREADS} \
      -b \
      -q ${MAPQ} \
      ${FILTER_FLAGS} \
      "${IN_BAM}" \
      -o filtered.bam
  else
    echo "POL2 | BAM | Filtering and sorting BAM (MAPQ≥${MAPQ})..."
    samtools view \
      -@ ${THREADS} \
      -b \
      -q ${MAPQ} \
      ${FILTER_FLAGS} \
      "${IN_BAM}" | \
    samtools sort -@ ${THREADS} -o filtered.bam
  fi

  BAM_END=$(date +%s)
  BAM_TIME=$((BAM_END - BAM_START))
  
  echo "POL2 | BAM | Indexing filtered BAM..."
  samtools index -@ ${THREADS} filtered.bam

  FILT_SIZE=$(stat -c%s filtered.bam 2>/dev/null || stat -f%z filtered.bam 2>/dev/null || echo "unknown")
  FILT_READS=$(samtools view -c filtered.bam)
  
  echo "POL2 | BAM | Filtered BAM: ${FILT_SIZE} bytes (${FILT_READS} reads)"
  echo "POL2 | BAM | Processing time: ${BAM_TIME}s"

  ###########################################################################
  # 5) CALCULATE GENE METRICS AND PAUSING INDEX
  ###########################################################################

  echo "POL2 | GENES | Calculating per-gene metrics..."
  echo "POL2 | GENES | This may take several minutes for large gene sets..."

  GENES_START=$(date +%s)

  # Build optional flags
  FAIL_FLAG=""
  [[ ${FAIL_IF_NO_GENES} -eq 1 ]] && FAIL_FLAG="--fail-if-empty true"

  set +e
  python3 "${CALC_SCRIPT}" \
    --bam filtered.bam \
    --gtf "${GTF_FILE}" \
    --tss-win ${TSS_WIN} \
    --body-offset-min ${BODY_OFFSET_MIN} \
    --body-offset-frac ${BODY_OFFSET_FRAC} \
    --feature-types "${FEATURE_TYPES}" \
    --out-pausing pausing_index.tsv \
    --out-genes pol2_gene_metrics.tsv \
    --out-qc pol2_qc.json \
    --threads ${THREADS} \
    ${FAIL_FLAG}
  
  GENES_RC=$?
  set -e

  GENES_END=$(date +%s)
  GENES_TIME=$((GENES_END - GENES_START))

  echo "POL2 | GENES | Processing completed in ${GENES_TIME}s"

  # Handle failures
  if [[ ${GENES_RC} -ne 0 ]]; then
    echo "POL2 | ERROR | Gene metrics calculation failed with exit code ${GENES_RC}"
    
    if [[ ${FAIL_IF_NO_GENES} -eq 1 ]]; then
      echo "POL2 | ERROR | fail_if_no_genes=true, exiting"
      exit 1
    else
      echo "POL2 | WARNING | Creating empty output files"
    fi
  fi

  # Ensure output files exist
  if [[ ! -s pol2_gene_metrics.tsv ]]; then
    cat > pol2_gene_metrics.tsv <<GENESEOF
gene_id gene_name chrom strand  tss_bp  tss_lo  tss_hi  tss_width body_lo body_hi body_len  tss_count tss_cpm tss_density_per_bp  body_count  body_cpm  body_density_per_bp pi_raw  pi_len_norm is_truncated
GENESEOF
  fi

  if [[ ! -s pausing_index.tsv ]]; then
    cat > pausing_index.tsv <<PAUSINGEOF
gene_id	chrom	strand	tss_count	gene_body_count	pi_raw	pi_len_norm	is_truncated
PAUSINGEOF
  fi

  # Count results
  GENE_COUNT=$(tail -n +2 pol2_gene_metrics.tsv 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  PAUSING_COUNT=$(tail -n +2 pausing_index.tsv 2>/dev/null | wc -l | tr -d ' ' || echo 0)

  echo "POL2 | GENES | Gene metrics: ${GENE_COUNT} genes"
  echo "POL2 | GENES | Pausing indices: ${PAUSING_COUNT} genes"

  ###########################################################################
  # 6) CREATE README
  ###########################################################################

  echo "POL2 | README | Creating documentation..."

  cat > README_pol2_metrics.txt <<'DOCEOF'
================================================================================
POL-II METRICS — !{sid}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  RNA Polymerase II metrics quantifying transcription at multiple levels:
  
  1. Density Metrics: Signal per functional region (promoter, gene body, etc.)
  2. Gene Metrics: Per-gene TSS and body coverage
  3. Pausing Index: Ratio indicating promoter-proximal pausing

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample:     !{sid}
  Condition:  !{cond}
  Timepoint:  !{tp}
  Replicate:  !{rep}

DENSITY METRICS (from normalized tracks)
────────────────────────────────────────────────────────────────────────────
  Method: Signal quantification from normalized 3' coverage
  
  Normalization: ${NORM_METHOD}
    • siCPM: Spike-in normalized (preferred when available)
    • CPM: Standard library size normalization (fallback)
  
  Processing:
    1. Take absolute values: |positive| and |negative| tracks
    2. Merge strands: combined signal = |pos| + |neg|
    3. Map to regions: bedtools map -o sum
    4. Quantify per functional region
  
  Output: pol2_density.tsv
    Columns: chr, start, end, name, signal, norm_method
  
  Results: $([ -s pol2_density.tsv ] && echo "$(tail -n +2 pol2_density.tsv | wc -l) regions" || echo "No data")

GENE METRICS (from BAM alignments)
────────────────────────────────────────────────────────────────────────────
  Method: Read counting in TSS and gene body windows
  
  BAM Filtering:
    • MAPQ threshold: ≥${MAPQ}
    • Duplicates: $([ ${DEDUP_ENABLED} -eq 1 ] && echo "Removed" || echo "Retained")
    • Unmapped reads: Excluded
    • Filtered reads: ${FILT_READS}
  
  TSS Window:
    • Definition: ±${TSS_WIN} bp around transcription start site
    • Purpose: Measures promoter occupancy
    • Total width: $((TSS_WIN * 2)) bp
  
  Gene Body:
    • Start: max(${BODY_OFFSET_MIN} bp, ${BODY_OFFSET_FRAC} × gene_length)
    • End: Gene end (TES or last exon)
    • Purpose: Measures productive elongation
    • Excludes promoter-proximal region
  
  Metrics Calculated:
    • TSS count: Reads in TSS window
    • TSS CPM: Normalized TSS count
    • TSS density: Reads per base pair in TSS
    • Body count: Reads in gene body
    • Body CPM: Normalized body count
    • Body density: Reads per base pair in body
  
  Output: pol2_gene_metrics.tsv
    Contains per-gene coverage and density values
  
  Results: ${GENE_COUNT} genes analyzed

PAUSING INDEX
────────────────────────────────────────────────────────────────────────────
  Definition:
    PI = (TSS_density / TSS_width) / (Body_density / Body_width)
    
    Or equivalently:
    PI = (TSS_count / TSS_width) / (Body_count / Body_length)
  
  Interpretation:
    • PI > 1.5:  Strong promoter-proximal pausing
                 Pol-II accumulates at TSS, limited elongation
                 Common in developmental/stress-response genes
    
    • PI ≈ 1.0:  Balanced distribution
                 Normal transcription dynamics
    
    • PI < 1.0:  Productive elongation
                 Efficient transition to elongation
                 Common in housekeeping genes
  
  Two Variants:
    • pi_raw: Direct ratio (more interpretable)
    • pi_len_norm: Length-normalized (accounts for window sizes)
  
  Truncation:
    Genes marked "is_truncated" if body region couldn't be
    properly defined (e.g., overlapping regions)
  
  Output: pausing_index.tsv
    Contains per-gene pausing indices
  
  Results: ${PAUSING_COUNT} genes with pausing indices

FILES
────────────────────────────────────────────────────────────────────────────
  pol2_gene_metrics.tsv     — Per-gene TSS and body metrics
  pausing_index.tsv         — Per-gene pausing indices
  pol2_density.tsv          — Signal per functional region
  pol2_qc.json              — Quality control summary (optional)
  README_pol2_metrics.txt   — This documentation
  pol2_metrics.log          — Processing log

FILE FORMATS
────────────────────────────────────────────────────────────────────────────

pol2_gene_metrics.tsv:
  gene_id              — Gene identifier
  gene_name            — Gene symbol/name
  chrom                — Chromosome
  strand               — Strand (+/-)
  tss_bp               — TSS coordinate
  tss_lo, tss_hi       — TSS window bounds
  tss_width            — TSS window width
  body_lo, body_hi     — Gene body bounds
  body_len             — Gene body length
  tss_count            — Reads in TSS
  tss_cpm              — TSS CPM
  tss_density_per_bp   — TSS reads per bp
  body_count           — Reads in body
  body_cpm             — Body CPM
  body_density_per_bp  — Body reads per bp
  pi_raw               — Raw pausing index
  pi_len_norm          — Length-normalized PI
  is_truncated         — Boolean flag

pausing_index.tsv:
  gene_id              — Gene identifier
  chrom                — Chromosome
  strand               — Strand
  tss_count            — TSS read count
  gene_body_count      — Body read count
  pausing_index        — Calculated PI
  is_truncated         — Boolean flag

pol2_density.tsv:
  chr                  — Chromosome
  start                — Region start
  end                  — Region end
  name                 — Region name/category
  signal               — Quantified signal
  norm_method          — Normalization (CPM or siCPM)

pol2_qc.json:
  JSON format with:
    • Total genes processed
    • Genes with data
    • Mean pausing index
    • Median pausing index
    • Coverage statistics

QUALITY CONTROL
────────────────────────────────────────────────────────────────────────────

Expected Values:
  • Pausing Index: 0.5 - 3.0 for most genes
  • TSS density: > 0 for active genes
  • Body density: > 0 for transcribed genes
  • Filtered reads: >1M for good coverage

Red Flags:
  • Very low gene count: Check GTF compatibility, feature types
  • All PI = 0 or NaN: Check BAM quality, TSS coordinates
  • Extreme PI values (>10): Check gene annotations
  • Most genes truncated: Body offset too large

Troubleshooting:
  • Zero genes: Check feature_types parameter, GTF format
  • Low coverage: Increase sequencing depth, check alignment
  • High truncation: Reduce body_offset_min or body_offset_frac
  • Missing density: Check functional regions availability

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These metrics are used for:
  
  1. Differential Analysis:
     - Compare pausing indices across conditions
     - Identify genes with changed Pol-II dynamics
  
  2. Quality Control:
     - Assess library quality
     - Validate expected expression patterns
  
  3. Mechanistic Studies:
     - Identify paused promoters
     - Study elongation regulation
     - Compare promoter vs. elongation changes
  
  4. Integration:
     - Correlate with ChIP-seq (Pol-II, NELF, DSIF)
     - Compare with nascent RNA-seq
     - Integrate with gene expression data

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  MAPQ threshold:      ${MAPQ}
  Remove duplicates:   $([ ${DEDUP_ENABLED} -eq 1 ] && echo "Yes" || echo "No")
  TSS window:          ±${TSS_WIN} bp
  Body offset min:     ${BODY_OFFSET_MIN} bp
  Body offset frac:    ${BODY_OFFSET_FRAC}
  Feature types:       ${FEATURE_TYPES}
  Normalization:       ${NORM_METHOD}
  
PROCESSING TIMES
────────────────────────────────────────────────────────────────────────────
  BAM filtering:       ${BAM_TIME}s
  Gene metrics:        ${GENES_TIME}s
  Total:               $((BAM_TIME + GENES_TIME))s

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • Uses bedtools coverage for efficient counting
  • Coordinate-sorted BAMs processed faster
  • Linear algebra threading limited to prevent oversubscription
  • Handles both standard GTF and GFF3 formats
  • Gene body offset prevents TSS/body overlap

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: !{sid}
  Module: 11_calculate_pol2_metrics

================================================================================
DOCEOF

  echo "POL2 | README | Documentation created"

  ###########################################################################
  # 7) VALIDATION
  ###########################################################################

  echo "POL2 | VALIDATE | Validating outputs..."

  VALIDATION_OK=1

  # Check required outputs
  for FILE in pol2_gene_metrics.tsv pausing_index.tsv pol2_density.tsv; do
    if [[ ! -e "${FILE}" ]]; then
      echo "POL2 | ERROR | Missing output file: ${FILE}"
      VALIDATION_OK=0
    fi
  done

  # Validate file formats
  if [[ -s pol2_gene_metrics.tsv ]]; then
    GENES_COLS=$(head -1 pol2_gene_metrics.tsv | awk -F'\t' '{print NF}')
    if [[ ${GENES_COLS} -ne 19 ]]; then
      echo "POL2 | WARNING | Gene metrics has ${GENES_COLS} columns, expected 19"
    fi
  fi

  if [[ -s pausing_index.tsv ]]; then
    PAUSING_COLS=$(head -1 pausing_index.tsv | awk -F'\t' '{print NF}')
    if [[ ${PAUSING_COLS} -ne 7 ]]; then
      echo "POL2 | WARNING | Pausing index has ${PAUSING_COLS} columns, expected 7"
    fi
  fi

  if [[ -s pol2_density.tsv ]]; then
    DENSITY_COLS=$(head -1 pol2_density.tsv | awk -F'\t' '{print NF}')
    if [[ ${DENSITY_COLS} -ne 6 ]]; then
      echo "POL2 | WARNING | Density has ${DENSITY_COLS} columns, expected 6"
    fi
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "POL2 | ERROR | Validation failed"
    exit 1
  fi

  echo "POL2 | VALIDATE | All outputs validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  # Calculate some statistics if possible
  if [[ ${GENE_COUNT} -gt 0 && -s pausing_index.tsv ]]; then
    # Calculate median PI (excluding header and truncated)
    # Use pi_len_norm (col 7) if 8-column format, else pi_raw (col 6); exclude truncated (last col)
    MEDIAN_PI=$(tail -n +2 pausing_index.tsv | \
                awk -F'\t' '$NF!="1" && $NF!="True" && ((NF>=8 && $7!="NA" && $7+0>0) || (NF==7 && $6!="NA" && $6+0>0)) {print (NF>=8?$7:$6)}' | \
                sort -n | \
                awk '{a[NR]=$1} END{print (NR%2==1)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}' || echo "NA")
  else
    MEDIAN_PI="NA"
  fi

  echo "────────────────────────────────────────────────────────────────────────"
  echo "POL2 | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "POL2 | SUMMARY | Normalization: ${NORM_METHOD}"
  echo "POL2 | SUMMARY | Filtered reads: ${FILT_READS}"
  echo "POL2 | SUMMARY | Genes analyzed: ${GENE_COUNT}"
  echo "POL2 | SUMMARY | Pausing indices: ${PAUSING_COUNT}"
  echo "POL2 | SUMMARY | Median PI: ${MEDIAN_PI}"
  echo "POL2 | SUMMARY | Density regions: $([ -s pol2_density.tsv ] && tail -n +2 pol2_density.tsv | wc -l || echo 0)"
  echo "POL2 | SUMMARY | Processing time: $((BAM_TIME + GENES_TIME))s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "POL2 | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}