// ============================================================================
// normalize_coverage_tracks.nf — Track Normalization (CPM and Spike-in CPM)
// ============================================================================
//
// Purpose:
//   Normalizes raw coverage tracks to CPM and spike-in CPM (siCPM)
//
// Features:
//   • CPM normalization: Counts per million mapped reads
//   • siCPM normalization: Spike-in normalized CPM
//   • Processes both main and allMap BAM tracks
//   • Handles 3' tracks (always) and 5' tracks (always)
//   • Single-pass scaling for efficiency
//   • Optional BigWig generation
//   • Comprehensive manifest of all outputs
//
// Normalization Methods:
//   CPM   = (count / sample_reads) × 1,000,000
//   siCPM = (count / sample_spike) × (control_spike / control_reads) × 1,000,000
//
// siCPM Control Selection:
//   1. Primary: condition matches params.control_label AND replicate == 1
//   2. Fallback: First sample with spike_reads > 0
//
// Inputs:
//   tuple(sample_id, pos3_bg, neg3_bg, pos5_bg, neg5_bg,
//         allmap3p_pos_bg, allmap3p_neg_bg, allmap5p_pos_bg, allmap5p_neg_bg,
//         condition, timepoint, replicate, counts_master_tsv, genes_unused)
//   path(genome_fa) : Genome FASTA (for chromosome sizes)
//
// Outputs:
//   ${params.output_dir}/05_normalized_tracks/${sample_id}/
//     ├── cpm/
//     │   ├── 3p/${sample_id}.3p.pos.cpm.bedgraph
//     │   ├── 3p/${sample_id}.3p.neg.cpm.bedgraph
//     │   ├── 3p/${sample_id}.3p.pos.cpm.bw
//     │   ├── 3p/${sample_id}.3p.neg.cpm.bw
//     │   ├── 3p/${sample_id}.allMap.3p.*.{bedgraph,bw} (if emit_allmap)
//     │   └── 5p/ (if emit_5p or auto-detected) [same structure]
//     ├── sicpm/ (if spike-in available)
//     │   ├── 3p/${sample_id}.3p.pos.sicpm.{bedgraph,bw}
//     │   └── 5p/ ...
//     ├── gene_end/ (if norm.gene_end_method != none)
//     │   ├── 3p/${sample_id}.3p.pos.gene_end.{bedgraph,bw}
//     │   └── 5p/ ...
//     ├── normalization_factors.tsv
//     ├── tracks_manifest.tsv
//     ├── README_normalization.txt
//     └── normalize_coverage_tracks.log
//
// Parameters:
//   params.norm.emit_bw            : Generate BigWig files (default: true)
//   params.norm.emit_sicpm         : Generate siCPM tracks (default: true)
//   params.norm.emit_allmap        : Include allMap tracks (default: true)
//   params.norm.emit_5p            : Force 5' track generation (default: auto)
//   params.norm.gene_end_method    : none | tes_window | gene_body (default: none)
//   params.norm.gene_end_window    : bp upstream of TES to sample (default: 500)
//   params.norm.gene_end_min_reads : min reads for gene body inclusion (default: 10)
//   params.force_sort_bedgraph     : Sort before BigWig (default: false)
//   params.control_label           : Control condition name (default: "CTRL")
//   params.norm.timeout_bw         : BigWig timeout seconds (default: 900)
//
// ============================================================================

nextflow.enable.dsl = 2

process normalize_coverage_tracks {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/05_normalized_tracks/${sample_id}",
             mode: params.publish_mode,
             overwrite: true,
             saveAs: { filename ->
               def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
               // Skip bedGraphs when output.bedgraph: false (BigWigs sufficient for genome browsers)
               if (params.get('output')?.get('bedgraph') == false && name.endsWith('.bedgraph')) return null
               return name
             }

  conda (params.conda_norm ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  // NOTE: bedGraph inputs MUST be path() not val() for proper Docker staging
  input:
    tuple val(sample_id),
          path(pos3_bg),
          path(neg3_bg),
          path(pos5_bg),
          path(neg5_bg),
          path(am3p_pos_bg),
          path(am3p_neg_bg),
          path(am5p_pos_bg),
          path(am5p_neg_bg),
          val(condition), val(timepoint), val(replicate),
          path(counts_master_tsv),
          val(genes_unused)
    path genome_fa
    path tes_bed      // TES BED file for gene-end normalization (optional; can be empty file)

  // ── Outputs ───────────────────────────────────────────────────────────────
  //
  // Output directory structure (per-method subfolders):
  //   05_normalized_tracks/{sample_id}/
  //     ├── cpm/3p/       — CPM-normalized 3' tracks
  //     ├── cpm/5p/       — CPM-normalized 5' tracks (if enabled)
  //     ├── sicpm/3p/     — siCPM-normalized 3' tracks (if spike-in available)
  //     ├── sicpm/5p/     — siCPM-normalized 5' tracks (if enabled + spike-in)
  //     └── gene_end/3p/  — Gene-end normalized 3' tracks (if method != none)
  //         gene_end/5p/  — Gene-end normalized 5' tracks (if enabled + method != none)
  //
  output:
    // Legacy compatibility tuple (paths updated to new cpm/ subfolder)
    tuple val(sample_id),
          path("cpm/3p/${sample_id}_pos3_cpm.bedgraph"),
          path("cpm/3p/${sample_id}_neg3_cpm.bedgraph"),
          path("normalization_factors.tsv"),
          val(condition), val(timepoint), val(replicate),
          emit: norm_tuple

    // Main CPM tracks
    path "cpm/3p/${sample_id}.3p.pos.cpm.bedgraph", emit: pos3_cpm_bg
    path "cpm/3p/${sample_id}.3p.neg.cpm.bedgraph", emit: neg3_cpm_bg

    // Documentation and manifest
    path "README_normalization.txt", emit: readme
    path "tracks_manifest.tsv",      emit: manifest
    path "normalization_factors.tsv"

    // Full directory trees (per-method)
    path "cpm/3p/**",                              emit: tree_cpm_3p
    path "cpm/5p/**",  optional: true,             emit: tree_cpm_5p
    path "sicpm/3p/**", optional: true,            emit: tree_sicpm_3p
    path "sicpm/5p/**", optional: true,            emit: tree_sicpm_5p
    path "gene_end/3p/**", optional: true,         emit: tree_gene_end_3p
    path "gene_end/5p/**", optional: true,         emit: tree_gene_end_5p

    // Log
    path "normalize_coverage_tracks.log", emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  normEmitBw      = params.norm?.emit_bw != null ? params.norm.emit_bw : true
  normEmitSicpm   = params.norm?.emit_sicpm != null ? params.norm.emit_sicpm : true
  normEmitAllmap  = params.norm?.emit_allmap != null ? params.norm.emit_allmap : true
  normTimeoutBw   = params.norm?.timeout_bw ?: 900
  normEmit5p      = params.norm?.emit_5p
  normForceSortBedgraph  = params.force_sort_bedgraph ? 'true' : 'false'
  normControlLabel       = (params.control_label ?: 'CTRL').toString()
  normGeneEndMethod      = (params.norm?.gene_end_method ?: 'none').toString()
  normGeneEndWindow      = params.norm?.gene_end_window    ?: 500
  normGeneEndMinReads    = params.norm?.gene_end_min_reads ?: 10
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a normalize_coverage_tracks.log)
  exec 2> >(tee -a normalize_coverage_tracks.log >&2)

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
  trap 'tracktx_error "normalize_coverage_tracks" "Unexpected process failure" "Check normalize_coverage_tracks.log in work dir"' ERR

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "NORMALIZE | START | sample=!{sample_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sample_id}"
  CONDITION="!{condition}"
  TIMEPOINT="!{timepoint}"
  REPLICATE="!{replicate}"
  THREADS=!{task.cpus}
  
  COUNTS_MASTER="!{counts_master_tsv}"
  GENOME_FA="!{genome_fa}"
  
  # Input bedGraphs
  POS3="!{pos3_bg}"
  NEG3="!{neg3_bg}"
  POS5="!{pos5_bg}"
  NEG5="!{neg5_bg}"
  AM3P_POS="!{am3p_pos_bg}"
  AM3P_NEG="!{am3p_neg_bg}"
  AM5P_POS="!{am5p_pos_bg}"
  AM5P_NEG="!{am5p_neg_bg}"

  # Feature flags
  EMIT_BW=$([[ "!{normEmitBw}" == "false" ]] && echo 0 || echo 1)
  EMIT_SICPM=$([[ "!{normEmitSicpm}" == "false" ]] && echo 0 || echo 1)
  EMIT_ALLMAP=$([[ "!{normEmitAllmap}" == "false" ]] && echo 0 || echo 1)
  FORCE_SORT=$([[ "!{normForceSortBedgraph}" == "true" ]] && echo 1 || echo 0)

  CONTROL_LABEL="!{normControlLabel}"
  TIMEOUT_BW=!{normTimeoutBw}
  GENE_END_METHOD="!{normGeneEndMethod}"
  GENE_END_WINDOW=!{normGeneEndWindow}
  GENE_END_MIN_READS=!{normGeneEndMinReads}
  TES_BED="!{tes_bed}"

  EMIT_GENE_END=0
  [[ "${GENE_END_METHOD}" != "none" && "${GENE_END_METHOD}" != "" ]] && EMIT_GENE_END=1

  # Auto-detect 5' track generation
  EMIT_5P_SETTING="!{normEmit5p}"
  if [[ "${EMIT_5P_SETTING}" == "true" ]]; then
    EMIT_5P=1
  elif [[ "${EMIT_5P_SETTING}" == "false" ]]; then
    EMIT_5P=0
  else
    # Auto mode: enable if 5' inputs exist
    if [[ -e "${POS5}" || -e "${NEG5}" ]]; then
      EMIT_5P=1
    else
      EMIT_5P=0
    fi
  fi

  echo "NORMALIZE | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "NORMALIZE | CONFIG | Condition: ${CONDITION}"
  echo "NORMALIZE | CONFIG | Timepoint: ${TIMEPOINT}"
  echo "NORMALIZE | CONFIG | Replicate: ${REPLICATE}"
  echo "NORMALIZE | CONFIG | Threads: ${THREADS}"
  echo "NORMALIZE | CONFIG | Control label: ${CONTROL_LABEL}"
  echo "NORMALIZE | CONFIG | Emit BigWig: $([ ${EMIT_BW} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | Emit siCPM: $([ ${EMIT_SICPM} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | Emit allMap: $([ ${EMIT_ALLMAP} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | Emit 5' tracks: $([ ${EMIT_5P} -eq 1 ] && echo "yes" || echo "no (auto)")"
  echo "NORMALIZE | CONFIG | Force bedGraph sort: $([ ${FORCE_SORT} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | BigWig timeout: ${TIMEOUT_BW}s"
  echo "NORMALIZE | CONFIG | Gene-end method: ${GENE_END_METHOD}"
  [[ ${EMIT_GENE_END} -eq 1 ]] && echo "NORMALIZE | CONFIG | Gene-end window: ${GENE_END_WINDOW} bp | min reads: ${GENE_END_MIN_READS}"

  # Create output directories (per-method subfolders)
  mkdir -p cpm/3p cpm/5p sicpm/3p sicpm/5p
  [[ ${EMIT_GENE_END} -eq 1 ]] && mkdir -p gene_end/3p gene_end/5p
  echo "NORMALIZE | CONFIG | Output directories created (cpm/, sicpm/$([ ${EMIT_GENE_END} -eq 1 ] && echo ", gene_end/"))"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "NORMALIZE | VALIDATE | Checking input files..."

  # Validate counts master file
  if [[ ! -s "${COUNTS_MASTER}" ]]; then
    tracktx_error "normalize_coverage_tracks" "Counts master file missing or empty: ${COUNTS_MASTER}" "Check quantify_reads_per_gene produced counts TSV"
  fi

  COUNTS_SIZE=$(stat -c%s "${COUNTS_MASTER}" 2>/dev/null || stat -f%z "${COUNTS_MASTER}" 2>/dev/null || echo "unknown")
  COUNTS_LINES=$(wc -l < "${COUNTS_MASTER}" | tr -d ' ')
  echo "NORMALIZE | VALIDATE | Counts master: ${COUNTS_SIZE} bytes (${COUNTS_LINES} lines)"

  # Validate input bedGraphs
  INPUT_COUNT=0
  for BG in "${POS3}" "${NEG3}" "${AM3P_POS}" "${AM3P_NEG}"; do
    if [[ -s "${BG}" ]]; then
      INPUT_COUNT=$((INPUT_COUNT + 1))
      BG_SIZE=$(stat -c%s "${BG}" 2>/dev/null || stat -f%z "${BG}" 2>/dev/null || echo "unknown")
      echo "NORMALIZE | VALIDATE | Input bedGraph: $(basename ${BG}) (${BG_SIZE} bytes)"
    fi
  done

  echo "NORMALIZE | VALIDATE | Found ${INPUT_COUNT} non-empty input bedGraphs"

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  else
    PYTHON_CMD="python3"
  fi

  # Validate tools
  if ! ${PYTHON_CMD} --version >/dev/null 2>&1; then
    tracktx_error "normalize_coverage_tracks" "Python not found (tried: ${PYTHON_CMD})" "Use -profile docker or install Python"
  fi
  echo "NORMALIZE | VALIDATE | python: $(${PYTHON_CMD} --version 2>&1)"
  if ! command -v awk >/dev/null 2>&1; then
    tracktx_error "normalize_coverage_tracks" "Required tool not found: awk" "Use -profile docker"
  fi
  echo "NORMALIZE | VALIDATE | awk: $(which awk)"
  if [[ ${EMIT_BW} -eq 1 ]] && ! command -v bedGraphToBigWig >/dev/null 2>&1; then
    tracktx_error "normalize_coverage_tracks" "bedGraphToBigWig not found (required for BigWig)" "Install UCSC tools or use -profile docker"
  fi
  if [[ ${EMIT_BW} -eq 1 ]]; then
    echo "NORMALIZE | VALIDATE | bedGraphToBigWig: $(which bedGraphToBigWig)"
  fi

  ###########################################################################
  # 3) COMPUTE NORMALIZATION FACTORS
  ###########################################################################

  echo "NORMALIZE | FACTORS | Computing CPM and siCPM scaling factors..."

  # Python script to compute factors from counts master
  ${PYTHON_CMD} - "${COUNTS_MASTER}" "${SAMPLE_ID}" "${CONTROL_LABEL}" > factors.tmp <<'PYSCRIPT'
import sys
import csv

counts_file = sys.argv[1]
sample_id = sys.argv[2]
control_label = sys.argv[3].strip().lower()

# Read counts file
with open(counts_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)

# Helper function for case-insensitive comparison
def normalize_str(s):
    return (s or '').strip().lower()

# Find current sample
sample_row = None
for row in rows:
    if row.get('sample') == sample_id:
        sample_row = row
        break

if sample_row is None:
    print("0.0000000000\t0.0000000000")
    sys.exit(0)

# Extract sample counts
sample_main = int(sample_row.get('main_reads', 0))
sample_spike = int(sample_row.get('spike_reads', 0))

# Compute CPM factor
if sample_main > 0:
    fac_cpm = 1_000_000.0 / sample_main
else:
    fac_cpm = 0.0

# Find control sample for siCPM
control_row = None

# First try: match control_label and replicate 1
for row in rows:
    cond = normalize_str(row.get('condition', ''))
    rep = str(row.get('replicate', '')).strip()
    if cond == control_label and rep in ('1', 'r1', 'R1'):
        control_row = row
        break

# Fallback: first sample with spike_reads > 0
if control_row is None:
    for row in rows:
        if int(row.get('spike_reads', 0)) > 0:
            control_row = row
            break

# Compute siCPM factor
if control_row and sample_spike > 0:
    control_main = int(control_row.get('main_reads', 0))
    control_spike = int(control_row.get('spike_reads', 0))
    
    if control_spike > 0 and control_main > 0:
        fac_sicpm = (control_spike / float(sample_spike)) * (1_000_000.0 / control_main)
    else:
        fac_sicpm = 0.0
else:
    fac_sicpm = 0.0

print(f"{fac_cpm:.10f}\t{fac_sicpm:.10f}")
PYSCRIPT

  # Read computed factors
  read -r FAC_CPM FAC_SICPM < factors.tmp

  echo "NORMALIZE | FACTORS | CPM factor: ${FAC_CPM}"
  echo "NORMALIZE | FACTORS | siCPM factor: ${FAC_SICPM}"

  # Validate CPM factor
  if awk -v x="${FAC_CPM}" 'BEGIN{exit (x>0?0:1)}'; then
    echo "NORMALIZE | FACTORS | CPM normalization enabled"
  else
    tracktx_error "normalize_coverage_tracks" "Cannot compute CPM (sample reads = 0)" "Check quantify_reads_per_gene output and counts TSV"
  fi

  # Check siCPM availability
  if awk -v x="${FAC_SICPM}" 'BEGIN{exit (x>0?0:1)}'; then
    echo "NORMALIZE | FACTORS | siCPM normalization enabled"
    SICPM_AVAILABLE=1
  else
    echo "NORMALIZE | FACTORS | WARNING: siCPM disabled (factor = 0)"
    echo "NORMALIZE | FACTORS | Possible causes: no spike-in reads, no control sample"
    SICPM_AVAILABLE=0
  fi

  ###########################################################################
  # 4) PREPARE GENOME SIZES
  ###########################################################################

  echo "NORMALIZE | GENOME | Preparing chromosome sizes..."

  # Check for existing .fai
  GENOME_FA_SRC="!{genome_fa.toString()}"
  if [[ -s "${GENOME_FA_SRC}.fai" && ! -s "${GENOME_FA}.fai" ]]; then
    ln -sf "${GENOME_FA_SRC}.fai" "${GENOME_FA}.fai" 2>/dev/null || true
  fi

  # Create genome.sizes
  if [[ -s "${GENOME_FA}.fai" ]]; then
    echo "NORMALIZE | GENOME | Using existing FASTA index"
    cut -f1,2 "${GENOME_FA}.fai" > genome.sizes
  else
    echo "NORMALIZE | GENOME | Creating FASTA index..."
    samtools faidx "${GENOME_FA}"
    cut -f1,2 "${GENOME_FA}.fai" > genome.sizes
  fi

  # Sort and validate genome.sizes
  LC_ALL=C sort -k1,1 -u -o genome.sizes genome.sizes

  CHR_COUNT=$(wc -l < genome.sizes | tr -d ' ')
  TOTAL_SIZE=$(awk '{sum+=$2} END{print sum}' genome.sizes)
  echo "NORMALIZE | GENOME | Chromosomes: ${CHR_COUNT}"
  echo "NORMALIZE | GENOME | Total size: ${TOTAL_SIZE} bp"

  ###########################################################################
  # 5) HELPER FUNCTIONS
  ###########################################################################

  # Convert bedGraph to BigWig
  make_bigwig() {
    local bedgraph="$1"
    local bigwig="$2"
    
    if [[ ! -s "${bedgraph}" ]]; then
      echo "NORMALIZE | BIGWIG | Empty input, creating empty BigWig: $(basename ${bigwig})"
      : > "${bigwig}"
      return 0
    fi
    
    if [[ ${EMIT_BW} -eq 0 ]]; then
      : > "${bigwig}"
      return 0
    fi
    
    echo "NORMALIZE | BIGWIG | Converting: $(basename ${bedgraph}) → $(basename ${bigwig})"
    
    # Optional sorting
    if [[ ${FORCE_SORT} -eq 1 ]]; then
      echo "NORMALIZE | BIGWIG | Sorting bedGraph..."
      LC_ALL=C sort -k1,1 -k2,2n -o "${bedgraph}" "${bedgraph}"
    fi
    
    # Count lines
    LINE_COUNT=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "${bedgraph}")
    echo "NORMALIZE | BIGWIG | bedGraph lines: ${LINE_COUNT}"
    
    # Convert with timeout
    if timeout "${TIMEOUT_BW}" bedGraphToBigWig "${bedgraph}" genome.sizes "${bigwig}"; then
      BW_SIZE=$(stat -c%s "${bigwig}" 2>/dev/null || stat -f%z "${bigwig}" 2>/dev/null || echo "unknown")
      echo "NORMALIZE | BIGWIG | Created: $(basename ${bigwig}) (${BW_SIZE} bytes)"
    else
      echo "NORMALIZE | BIGWIG | WARNING: Conversion failed or timed out, creating empty BigWig"
      : > "${bigwig}"
    fi
  }

  # Normalize bedGraph (single-pass CPM + siCPM) — writes to per-method subfolders
  normalize_bedgraph() {
    local input_bg="$1"
    local end_label="$2"      # "3p" or "5p"
    local strand="$3"          # "pos" or "neg"
    local set_label="$4"       # "main" or "allMap"

    # Build output file names (method subfolder / end subfolder / filename)
    local set_suffix=""
    [[ "${set_label}" == "allMap" ]] && set_suffix=".allMap"

    local prefix="${SAMPLE_ID}${set_suffix}.${end_label}.${strand}"
    local out_cpm_bg="cpm/${end_label}/${prefix}.cpm.bedgraph"
    local out_cpm_bw="cpm/${end_label}/${prefix}.cpm.bw"
    local out_sicpm_bg="sicpm/${end_label}/${prefix}.sicpm.bedgraph"
    local out_sicpm_bw="sicpm/${end_label}/${prefix}.sicpm.bw"

    echo "NORMALIZE | SCALE | Processing: ${end_label} ${set_label} ${strand}"

    if [[ ! -s "${input_bg}" ]]; then
      echo "NORMALIZE | SCALE | WARNING: Input empty, creating empty outputs"
      : > "${out_cpm_bg}"
      : > "${out_cpm_bw}"
      : > "${out_sicpm_bg}"
      : > "${out_sicpm_bw}"
      return 0
    fi

    # Count input lines
    INPUT_LINES=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "${input_bg}")
    echo "NORMALIZE | SCALE | Input lines: ${INPUT_LINES}"

    # Normalize with awk (single pass for both CPM and siCPM)
    START_TIME=$(date +%s)

    if [[ ${EMIT_SICPM} -eq 1 && ${SICPM_AVAILABLE} -eq 1 ]]; then
      echo "NORMALIZE | SCALE | Writing CPM and siCPM..."
      awk -v fc="${FAC_CPM}" -v fs="${FAC_SICPM}" -v OFS='\t' \
          -v cpm_out="${out_cpm_bg}" -v sicpm_out="${out_sicpm_bg}" '
        (NF>=4) && ($0!~/^(track|browser|#)/) {
          print $1, $2, $3, $4 * fc   > cpm_out
          print $1, $2, $3, $4 * fs   > sicpm_out
        }
      ' "${input_bg}"
    else
      echo "NORMALIZE | SCALE | Writing CPM only..."
      awk -v fc="${FAC_CPM}" -v OFS='\t' '
        (NF>=4) && ($0!~/^(track|browser|#)/) { print $1, $2, $3, $4 * fc }
      ' "${input_bg}" > "${out_cpm_bg}"
      : > "${out_sicpm_bg}"
    fi

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    [[ ${ELAPSED} -lt 1 ]] && ELAPSED=1
    RATE=$((INPUT_LINES / ELAPSED))
    echo "NORMALIZE | SCALE | Normalization complete in ${ELAPSED}s (~${RATE} lines/s)"

    # Convert to BigWig
    make_bigwig "${out_cpm_bg}" "${out_cpm_bw}"
    if [[ -s "${out_sicpm_bg}" ]]; then
      make_bigwig "${out_sicpm_bg}" "${out_sicpm_bw}"
    else
      : > "${out_sicpm_bw}"
    fi

    # Add to manifest
    echo -e "${SAMPLE_ID}\t${end_label}\t${set_label}\t${strand}\tcpm\t${out_cpm_bg}" >> tracks_manifest.tsv
    [[ -s "${out_sicpm_bg}" ]] && \
      echo -e "${SAMPLE_ID}\t${end_label}\t${set_label}\t${strand}\tsicpm\t${out_sicpm_bg}" >> tracks_manifest.tsv

    echo "NORMALIZE | SCALE | Complete: ${end_label} ${set_label} ${strand}"
  }

  ###########################################################################
  # 5b) GENE-END NORMALIZATION HELPERS
  ###########################################################################
  #
  # Two methods:
  #   tes_window : Count signal in a window upstream of each TES.
  #                Window = [TES - GENE_END_WINDOW, TES - 100] bp (strand-aware).
  #                Scale factor = TARGET_SUM / sample_window_signal.
  #
  #   gene_body  : Compute median signal density in gene bodies
  #                (TSS+500 to TES-500 for genes with ≥ GENE_END_MIN_READS reads).
  #                Scale factor = TARGET_MEDIAN / sample_median.
  #
  # Both methods normalise the same raw bedGraphs; they only differ in how
  # the per-sample scaling scalar is derived.

  TARGET_GENE_END_SUM=1000000   # Target total signal (analogous to 1M in CPM)
  FAC_GENE_END=0                # Set below if method != none

  compute_gene_end_factor() {
    local input_bg="$1"   # The primary 3' pos bedGraph (strand-aware signal)
    echo "NORMALIZE | GENE_END | Computing gene-end factor (method: ${GENE_END_METHOD})..."

    if [[ ! -s "${TES_BED}" ]]; then
      echo "NORMALIZE | GENE_END | WARNING: TES BED file empty or missing — skipping gene-end normalization"
      FAC_GENE_END=0
      return 0
    fi
    if [[ ! -s "${input_bg}" ]]; then
      echo "NORMALIZE | GENE_END | WARNING: Input bedGraph empty — skipping gene-end normalization"
      FAC_GENE_END=0
      return 0
    fi

    if [[ "${GENE_END_METHOD}" == "tes_window" ]]; then
      # ── TES-window method ─────────────────────────────────────────────────
      # Build a BED of windows upstream of each TES (strand-aware).
      # Format of TES BED: chr start end gene_id score strand
      # For + strand: window = [TES_pos - GENE_END_WINDOW, TES_pos - 100]
      # For - strand: window = [TES_pos + 100, TES_pos + GENE_END_WINDOW]
      awk -v win="${GENE_END_WINDOW}" 'BEGIN{OFS="\t"} NF>=6 {
        chr=$1; start=$2; end=$3; name=$4; score=$5; strand=$6
        # TES position
        if (strand == "+") {
          tes = end
          w_start = tes - win
          w_end   = tes - 100
        } else {
          tes = start
          w_start = tes + 100
          w_end   = tes + win
        }
        if (w_start < 0) w_start = 0
        if (w_end > w_start + 50) print chr, w_start, w_end, name, score, strand
      }' "${TES_BED}" > tes_windows.bed

      WIN_COUNT=$(wc -l < tes_windows.bed | tr -d ' ')
      echo "NORMALIZE | GENE_END | TES windows: ${WIN_COUNT}"

      if [[ "${WIN_COUNT}" -lt 10 ]]; then
        echo "NORMALIZE | GENE_END | WARNING: Too few TES windows (< 10) — skipping"
        FAC_GENE_END=0
        return 0
      fi

      # Compute sum of bedGraph signal overlapping TES windows
      WINDOW_SUM=$(bedtools intersect \
        -a tes_windows.bed \
        -b "${input_bg}" \
        -wb 2>/dev/null \
      | awk '{
          # last field is the bedGraph score; weight by overlap width
          w = $3 - $2
          b_start = $NF-3; b_end = $(NF-2)
          overlap_start = ($2 > b_start) ? $2 : b_start
          overlap_end   = ($3 < b_end)   ? $3 : b_end
          ov = (overlap_end - overlap_start > 0) ? (overlap_end - overlap_start) : 0
          sum += $NF * ov
        } END { print (sum > 0 ? sum : 0) }')

      echo "NORMALIZE | GENE_END | TES window signal sum: ${WINDOW_SUM}"
      FAC_GENE_END=$(awk -v s="${WINDOW_SUM}" -v t="${TARGET_GENE_END_SUM}" \
        'BEGIN { print (s > 0) ? t/s : 0 }')

    elif [[ "${GENE_END_METHOD}" == "gene_body" ]]; then
      # ── Gene-body method ──────────────────────────────────────────────────
      # Use TES BED combined with the genes BED (approximated: require
      # gene length info). Build gene body windows: for each gene in TES BED,
      # extend [TSS+500, TES-500] using the feature coordinates in TES BED.
      # TES BED col 2/3 gives the gene span; skip short genes (<2000 bp).
      awk -v mr="${GENE_END_MIN_READS}" 'BEGIN{OFS="\t"} NF>=6 {
        chr=$1; gstart=$2; gend=$3; name=$4; score=$5; strand=$6
        glen = gend - gstart
        if (glen < 2000) next    # skip short genes
        if (strand == "+") {
          b_start = gstart + 500
          b_end   = gend   - 500
        } else {
          b_start = gstart + 500
          b_end   = gend   - 500
        }
        if (b_end > b_start + 100) print chr, b_start, b_end, name, score, strand
      }' "${TES_BED}" > gene_bodies.bed

      GENE_COUNT=$(wc -l < gene_bodies.bed | tr -d ' ')
      echo "NORMALIZE | GENE_END | Gene body windows: ${GENE_COUNT}"

      if [[ "${GENE_COUNT}" -lt 10 ]]; then
        echo "NORMALIZE | GENE_END | WARNING: Too few gene bodies (< 10) — skipping"
        FAC_GENE_END=0
        return 0
      fi

      # Compute per-gene body density (signal / bp) and take median
      MEDIAN_DENSITY=$(bedtools intersect \
        -a gene_bodies.bed \
        -b "${input_bg}" \
        -wb 2>/dev/null \
      | awk '{
          win_len = $3 - $2
          if (win_len < 1) next
          gene_key = $4
          # accumulate weighted signal per gene
          sig[gene_key]  += $NF * ($NF-2 > $2 ? ($NF-2 < $3 ? $NF-2 : $3) - ($NF-3 > $2 ? $NF-3 : $2) : 0)
          len[gene_key]   = win_len
        } END {
          n = 0
          for (g in sig) { dens[n++] = sig[g] / len[g] }
          # simple median
          for (i=0; i<n-1; i++) for (j=i+1; j<n; j++) if (dens[i]>dens[j]) {t=dens[i];dens[i]=dens[j];dens[j]=t}
          if (n > 0) print dens[int(n/2)]
          else print 0
        }')

      echo "NORMALIZE | GENE_END | Median gene body density: ${MEDIAN_DENSITY}"
      FAC_GENE_END=$(awk -v d="${MEDIAN_DENSITY}" -v t="${TARGET_GENE_END_SUM}" \
        'BEGIN { print (d > 0) ? t/d : 0 }')

    else
      FAC_GENE_END=0
    fi

    echo "NORMALIZE | GENE_END | Scale factor: ${FAC_GENE_END}"
  }

  # Apply gene-end normalization factor to a bedGraph → gene_end/{end_label}/
  normalize_bedgraph_gene_end() {
    local input_bg="$1"
    local end_label="$2"     # "3p" or "5p"
    local strand="$3"         # "pos" or "neg"
    local set_label="$4"      # "main" or "allMap"

    local set_suffix=""
    [[ "${set_label}" == "allMap" ]] && set_suffix=".allMap"

    local prefix="${SAMPLE_ID}${set_suffix}.${end_label}.${strand}"
    local out_bg="gene_end/${end_label}/${prefix}.gene_end.bedgraph"
    local out_bw="gene_end/${end_label}/${prefix}.gene_end.bw"

    echo "NORMALIZE | GENE_END | Applying to: ${end_label} ${set_label} ${strand}"

    if [[ ! -s "${input_bg}" || "${FAC_GENE_END}" == "0" ]]; then
      echo "NORMALIZE | GENE_END | WARNING: Skipping (empty input or zero factor)"
      : > "${out_bg}"
      : > "${out_bw}"
      return 0
    fi

    awk -v fg="${FAC_GENE_END}" -v OFS='\t' '
      (NF>=4) && ($0!~/^(track|browser|#)/) { print $1, $2, $3, $4 * fg }
    ' "${input_bg}" > "${out_bg}"

    make_bigwig "${out_bg}" "${out_bw}"
    echo -e "${SAMPLE_ID}\t${end_label}\t${set_label}\t${strand}\tgene_end\t${out_bg}" >> tracks_manifest.tsv
    echo "NORMALIZE | GENE_END | Complete: ${end_label} ${set_label} ${strand}"
  }

  ###########################################################################
  # 6) COMPUTE GENE-END NORMALIZATION FACTOR (if enabled)
  ###########################################################################

  if [[ ${EMIT_GENE_END} -eq 1 ]]; then
    echo "NORMALIZE | GENE_END | Computing gene-end scaling factor from positive 3' bedGraph..."
    compute_gene_end_factor "${POS3}"
  fi

  ###########################################################################
  # 7) NORMALIZE MAIN 3' TRACKS
  # pos and neg are independent — run in parallel
  ###########################################################################

  echo "NORMALIZE | MAIN3P | Normalizing main 3' tracks (pos + neg in parallel)..."

  # Initialize manifest before launching parallel jobs
  : > tracks_manifest.tsv

  pids=()
  normalize_bedgraph "${POS3}" "3p" "pos" "main" & pids+=($!)
  normalize_bedgraph "${NEG3}" "3p" "neg" "main" & pids+=($!)
  NORM_FAILED=0
  for pid in "${pids[@]}"; do wait "${pid}" || NORM_FAILED=1; done
  [[ ${NORM_FAILED} -ne 0 ]] && tracktx_error "normalize_coverage_tracks" "One or more main 3' normalization jobs failed" "Check normalize_coverage_tracks.log in work dir"

  # Gene-end: main 3' tracks
  if [[ ${EMIT_GENE_END} -eq 1 ]]; then
    mkdir -p gene_end/3p
    pids=()
    normalize_bedgraph_gene_end "${POS3}" "3p" "pos" "main" & pids+=($!)
    normalize_bedgraph_gene_end "${NEG3}" "3p" "neg" "main" & pids+=($!)
    for pid in "${pids[@]}"; do wait "${pid}" || true; done
  fi

  echo "NORMALIZE | MAIN3P | Main 3' tracks complete"

  ###########################################################################
  # 8) NORMALIZE MAIN 5' TRACKS (if enabled)
  # pos and neg are independent — run in parallel
  ###########################################################################

  if [[ ${EMIT_5P} -eq 1 ]]; then
    echo "NORMALIZE | MAIN5P | Normalizing main 5' tracks (pos + neg in parallel)..."

    pids=()
    normalize_bedgraph "${POS5}" "5p" "pos" "main" & pids+=($!)
    normalize_bedgraph "${NEG5}" "5p" "neg" "main" & pids+=($!)
    NORM_FAILED=0
    for pid in "${pids[@]}"; do wait "${pid}" || NORM_FAILED=1; done
    [[ ${NORM_FAILED} -ne 0 ]] && tracktx_error "normalize_coverage_tracks" "One or more main 5' normalization jobs failed" "Check normalize_coverage_tracks.log in work dir"

    # Gene-end: main 5' tracks
    if [[ ${EMIT_GENE_END} -eq 1 ]]; then
      mkdir -p gene_end/5p
      pids=()
      normalize_bedgraph_gene_end "${POS5}" "5p" "pos" "main" & pids+=($!)
      normalize_bedgraph_gene_end "${NEG5}" "5p" "neg" "main" & pids+=($!)
      for pid in "${pids[@]}"; do wait "${pid}" || true; done
    fi

    echo "NORMALIZE | MAIN5P | Main 5' tracks complete"
  else
    echo "NORMALIZE | MAIN5P | Skipping 5' tracks (not enabled)"
  fi

  ###########################################################################
  # 9) NORMALIZE ALLMAP 3' TRACKS (if enabled)
  # pos and neg are independent — run in parallel
  ###########################################################################

  if [[ ${EMIT_ALLMAP} -eq 1 ]]; then
    echo "NORMALIZE | ALLMAP3P | Normalizing allMap 3' tracks (pos + neg in parallel)..."

    pids=()
    normalize_bedgraph "${AM3P_POS}" "3p" "pos" "allMap" & pids+=($!)
    normalize_bedgraph "${AM3P_NEG}" "3p" "neg" "allMap" & pids+=($!)
    NORM_FAILED=0
    for pid in "${pids[@]}"; do wait "${pid}" || NORM_FAILED=1; done
    [[ ${NORM_FAILED} -ne 0 ]] && tracktx_error "normalize_coverage_tracks" "One or more allMap 3' normalization jobs failed" "Check normalize_coverage_tracks.log in work dir"

    # Gene-end: allMap 3' tracks
    if [[ ${EMIT_GENE_END} -eq 1 ]]; then
      pids=()
      normalize_bedgraph_gene_end "${AM3P_POS}" "3p" "pos" "allMap" & pids+=($!)
      normalize_bedgraph_gene_end "${AM3P_NEG}" "3p" "neg" "allMap" & pids+=($!)
      for pid in "${pids[@]}"; do wait "${pid}" || true; done
    fi

    echo "NORMALIZE | ALLMAP3P | AllMap 3' tracks complete"
  else
    echo "NORMALIZE | ALLMAP3P | Skipping allMap 3' tracks (not enabled)"
  fi

  ###########################################################################
  # 10) NORMALIZE ALLMAP 5' TRACKS (if both enabled)
  # pos and neg are independent — run in parallel
  ###########################################################################

  if [[ ${EMIT_ALLMAP} -eq 1 && ${EMIT_5P} -eq 1 ]]; then
    echo "NORMALIZE | ALLMAP5P | Normalizing allMap 5' tracks (pos + neg in parallel)..."

    pids=()
    normalize_bedgraph "${AM5P_POS}" "5p" "pos" "allMap" & pids+=($!)
    normalize_bedgraph "${AM5P_NEG}" "5p" "neg" "allMap" & pids+=($!)
    NORM_FAILED=0
    for pid in "${pids[@]}"; do wait "${pid}" || NORM_FAILED=1; done
    [[ ${NORM_FAILED} -ne 0 ]] && tracktx_error "normalize_coverage_tracks" "One or more allMap 5' normalization jobs failed" "Check normalize_coverage_tracks.log in work dir"

    # Gene-end: allMap 5' tracks
    if [[ ${EMIT_GENE_END} -eq 1 ]]; then
      pids=()
      normalize_bedgraph_gene_end "${AM5P_POS}" "5p" "pos" "allMap" & pids+=($!)
      normalize_bedgraph_gene_end "${AM5P_NEG}" "5p" "neg" "allMap" & pids+=($!)
      for pid in "${pids[@]}"; do wait "${pid}" || true; done
    fi

    echo "NORMALIZE | ALLMAP5P | AllMap 5' tracks complete"
  else
    echo "NORMALIZE | ALLMAP5P | Skipping allMap 5' tracks"
  fi

  ###########################################################################
  # 11) CREATE NORMALIZATION FACTORS FILE
  ###########################################################################

  echo "NORMALIZE | OUTPUT | Writing normalization factors..."

  {
    echo -e "method\tfactor"
    echo -e "CPM\t${FAC_CPM}"
    echo -e "siCPM\t${FAC_SICPM}"
    [[ ${EMIT_GENE_END} -eq 1 ]] && echo -e "gene_end_${GENE_END_METHOD}\t${FAC_GENE_END}"
  } > normalization_factors.tsv

  echo "NORMALIZE | OUTPUT | Normalization factors written"

  ###########################################################################
  # 12) CREATE LEGACY COPIES (updated to new cpm/ subfolder paths)
  ###########################################################################

  echo "NORMALIZE | LEGACY | Creating legacy compatibility copies..."

  # Legacy copies for downstream compatibility (symlinks not supported on NFS)
  cp -f "cpm/3p/${SAMPLE_ID}.3p.pos.cpm.bedgraph" "cpm/3p/${SAMPLE_ID}_pos3_cpm.bedgraph"
  cp -f "cpm/3p/${SAMPLE_ID}.3p.neg.cpm.bedgraph" "cpm/3p/${SAMPLE_ID}_neg3_cpm.bedgraph"

  echo "NORMALIZE | LEGACY | Legacy copies created"

  ###########################################################################
  # 12) CREATE README
  ###########################################################################

  echo "NORMALIZE | README | Creating documentation..."

  cat > README_normalization.txt <<'DOCEOF'
================================================================================
NORMALIZED TRACKS — !{sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Normalized coverage tracks for quantitative comparison across samples.
  
  Normalization methods:
    1. CPM (Counts Per Million): Standard library size normalization
    2. siCPM (Spike-in CPM): Accounts for global changes via spike-in RNA
    3. Gene-end (optional): TES-proximal or gene-body signal normalization

NORMALIZATION FORMULAS
────────────────────────────────────────────────────────────────────────────

CPM (Counts Per Million):
  CPM = (raw_count / sample_reads) × 1,000,000
  
  Use when:
    • Comparing samples with similar RNA content
    • Standard differential expression analysis
    • No spike-in RNA added

siCPM (Spike-in Normalized CPM):
  siCPM = (raw_count / sample_spike) × (control_spike / control_reads) × 1,000,000
  
  Use when:
    • Global transcription changes expected
    • Spike-in RNA added during library prep
    • Treatment affects overall RNA levels

CONTROL SELECTION (for siCPM)
────────────────────────────────────────────────────────────────────────────
  Priority 1: Condition = "${CONTROL_LABEL}" AND Replicate = 1
  Priority 2: First sample with spike_reads > 0 (fallback)

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample:     !{sample_id}
  Condition:  !{condition}
  Timepoint:  !{timepoint}
  Replicate:  !{replicate}

NORMALIZATION FACTORS
────────────────────────────────────────────────────────────────────────────
  CPM factor:   ${FAC_CPM}
  siCPM factor: ${FAC_SICPM}
  
  $([ ${SICPM_AVAILABLE} -eq 0 ] && echo "  Note: siCPM disabled (factor = 0)" || echo "  siCPM enabled and available")

FILES (per-method subfolder structure)
────────────────────────────────────────────────────────────────────────────

CPM Tracks (Always Generated):
  cpm/3p/!{sample_id}.3p.pos.cpm.bedgraph   — Positive strand CPM
  cpm/3p/!{sample_id}.3p.neg.cpm.bedgraph   — Negative strand CPM
  cpm/3p/!{sample_id}.3p.pos.cpm.bw         — BigWig format
  cpm/3p/!{sample_id}.3p.neg.cpm.bw         — BigWig format
  cpm/5p/ ...                                 — 5' end CPM (if enabled)

siCPM Tracks (if spike-in available):
  sicpm/3p/!{sample_id}.3p.pos.sicpm.bedgraph
  sicpm/3p/!{sample_id}.3p.neg.sicpm.bedgraph
  sicpm/3p/!{sample_id}.3p.pos.sicpm.bw
  sicpm/3p/!{sample_id}.3p.neg.sicpm.bw
  sicpm/5p/ ...                               — 5' end siCPM (if enabled)

Gene-End Tracks (if gene_end_method != none):
  gene_end/3p/!{sample_id}.3p.pos.gene_end.bedgraph
  gene_end/3p/!{sample_id}.3p.neg.gene_end.bedgraph
  gene_end/3p/!{sample_id}.3p.pos.gene_end.bw
  gene_end/5p/ ...

AllMap Tracks (if emit_allmap=true):
  cpm/3p/!{sample_id}.allMap.3p.*.cpm.bedgraph
  sicpm/3p/!{sample_id}.allMap.3p.*.sicpm.bedgraph
  gene_end/3p/!{sample_id}.allMap.3p.*.gene_end.bedgraph

Legacy Copies (for compatibility):
  cpm/3p/!{sample_id}_pos3_cpm.bedgraph → cpm/3p/!{sample_id}.3p.pos.cpm.bedgraph
  cpm/3p/!{sample_id}_neg3_cpm.bedgraph → cpm/3p/!{sample_id}.3p.neg.cpm.bedgraph

Metadata:
  normalization_factors.tsv — CPM and siCPM scaling factors
  tracks_manifest.tsv       — Complete list of all generated tracks
  README_normalization.txt  — This documentation
  normalize_coverage_tracks.log      — Processing log

FILE FORMAT
────────────────────────────────────────────────────────────────────────────
  bedGraph: chr<TAB>start<TAB>end<TAB>normalized_coverage
  BigWig:   Binary indexed format for genome browsers

PROCESSING NOTES
────────────────────────────────────────────────────────────────────────────
  • Single-pass scaling: CPM and siCPM computed together for efficiency
  • No coordinate clipping needed (validated by generate_coverage_tracks module)
  • Negative strand values preserved from upstream mirroring
  • BigWig timeout: ${TIMEOUT_BW} seconds
  • Optional bedGraph sorting: $([ ${FORCE_SORT} -eq 1 ] && echo "enabled" || echo "disabled")

TRACKS MANIFEST
────────────────────────────────────────────────────────────────────────────
  The tracks_manifest.tsv file contains:
    Column 1: sample       — Sample identifier
    Column 2: end          — 3p or 5p
    Column 3: set          — main or allMap
    Column 4: strand       — pos or neg
    Column 5: scale        — cpm or sicpm
    Column 6: path         — Relative path to bedGraph file

USAGE RECOMMENDATIONS
────────────────────────────────────────────────────────────────────────────

For Genome Browsers:
  • Use BigWig (.bw) files for visualization
  • Load both positive and negative strand tracks
  • Negative values indicate reverse strand coverage

For Downstream Analysis:
  • Use CPM tracks for standard comparisons
  • Use siCPM tracks when spike-in added
  • bedGraph format recommended for computational analysis
  • Manifest file helps programmatic access

For Differential Analysis:
  • CPM tracks suitable for DESeq2, edgeR input
  • siCPM accounts for global expression changes
  • Consider biological replicates for statistics

QUALITY CHECKS
────────────────────────────────────────────────────────────────────────────

Expected Behavior:
  • CPM values typically 0.01 - 1000 for expressed regions
  • siCPM values should be similar to CPM if no global changes
  • All tracks should have same number of intervals as input

Troubleshooting:
  • Empty outputs: Check input bedGraphs exist and non-empty
  • Zero siCPM factor: No spike-in reads or no control sample
  • BigWig conversion failure: Try --force_sort_bedgraph

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  Emit BigWig:      $([ ${EMIT_BW} -eq 1 ] && echo "Yes" || echo "No")
  Emit siCPM:       $([ ${EMIT_SICPM} -eq 1 ] && echo "Yes" || echo "No")
  Emit allMap:      $([ ${EMIT_ALLMAP} -eq 1 ] && echo "Yes" || echo "No")
  Emit 5' tracks:   $([ ${EMIT_5P} -eq 1 ] && echo "Yes" || echo "No (auto)")
  Control label:    ${CONTROL_LABEL}
  BigWig timeout:   ${TIMEOUT_BW}s
  Force sort:       $([ ${FORCE_SORT} -eq 1 ] && echo "Yes" || echo "No")

DOWNSTREAM MODULES
────────────────────────────────────────────────────────────────────────────
  These normalized tracks are used by:
  • Divergent transcription detection (allMap 3' tracks)
  • Functional region calling (main 3' tracks)
  • Pol-II metrics calculation (CPM or siCPM depending on spike-in)
  • Quality control and visualization

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: !{sample_id}
  Module: 08_normalize_coverage_tracks

================================================================================
DOCEOF

  echo "NORMALIZE | README | Documentation created"

  ###########################################################################
  # 13) VALIDATION AND SUMMARY
  ###########################################################################

  echo "NORMALIZE | VALIDATE | Validating outputs..."

  # Count output files (scan all method subfolders).
  # Each find is guarded with || true so a missing directory (e.g. gene_end/
  # when GENE_END_METHOD=none, or sicpm/ when no spike-in) does not propagate
  # a non-zero exit through the pipefail pipeline and trigger the ERR trap.
  CPM_BG_COUNT=$( (find cpm      -name "*.cpm.bedgraph"      -type f 2>/dev/null || true) | wc -l | tr -d ' ')
  SICPM_BG_COUNT=$((find sicpm   -name "*.sicpm.bedgraph"    -type f 2>/dev/null || true) | wc -l | tr -d ' ')
  CPM_BW_COUNT=$( (find cpm      -name "*.cpm.bw"            -type f 2>/dev/null || true) | wc -l | tr -d ' ')
  SICPM_BW_COUNT=$((find sicpm   -name "*.sicpm.bw"          -type f 2>/dev/null || true) | wc -l | tr -d ' ')
  GENE_END_BG_COUNT=$((find gene_end -name "*.gene_end.bedgraph" -type f 2>/dev/null || true) | wc -l | tr -d ' ')
  MANIFEST_LINES=$(wc -l < tracks_manifest.tsv | tr -d ' ')

  echo "NORMALIZE | VALIDATE | CPM bedGraphs: ${CPM_BG_COUNT}"
  echo "NORMALIZE | VALIDATE | siCPM bedGraphs: ${SICPM_BG_COUNT}"
  echo "NORMALIZE | VALIDATE | Gene-end bedGraphs: ${GENE_END_BG_COUNT}"
  echo "NORMALIZE | VALIDATE | CPM BigWigs: ${CPM_BW_COUNT}"
  echo "NORMALIZE | VALIDATE | siCPM BigWigs: ${SICPM_BW_COUNT}"
  echo "NORMALIZE | VALIDATE | Manifest entries: ${MANIFEST_LINES}"

  # Check critical files exist (updated paths)
  for file in \
    "cpm/3p/${SAMPLE_ID}.3p.pos.cpm.bedgraph" \
    "cpm/3p/${SAMPLE_ID}.3p.neg.cpm.bedgraph" \
    "normalization_factors.tsv" \
    "tracks_manifest.tsv"; do

    if [[ ! -s "${file}" ]]; then
      tracktx_error "normalize_coverage_tracks" "Missing or empty critical file: ${file}" "Check normalize_coverage_tracks.log in work dir"
    fi
  done

  echo "NORMALIZE | VALIDATE | All critical files present"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  # Calculate total output size
  TOTAL_SIZE=$(du -sh . 2>/dev/null | cut -f1 || echo "unknown")

  echo "────────────────────────────────────────────────────────────────────────"
  echo "NORMALIZE | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "NORMALIZE | SUMMARY | CPM factor: ${FAC_CPM}"
  echo "NORMALIZE | SUMMARY | siCPM factor: ${FAC_SICPM}"
  echo "NORMALIZE | SUMMARY | CPM bedGraphs: ${CPM_BG_COUNT}"
  echo "NORMALIZE | SUMMARY | siCPM bedGraphs: ${SICPM_BG_COUNT}"
  echo "NORMALIZE | SUMMARY | BigWig files: $((CPM_BW_COUNT + SICPM_BW_COUNT))"
  echo "NORMALIZE | SUMMARY | Total output: ${TOTAL_SIZE}"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "NORMALIZE | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}