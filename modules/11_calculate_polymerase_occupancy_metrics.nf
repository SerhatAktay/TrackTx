// ============================================================================
// calculate_polymerase_occupancy_metrics.nf — Polymerase Occupancy Metrics Calculation
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
//   ${params.output_dir}/08_pol_metrics/${sample_id}/
//     ├── pol_gene_metrics.tsv    — Per-gene TSS and body metrics
//     ├── pausing_index.tsv        — Per-gene pausing indices
//     ├── pol_density.tsv         — Signal per functional region
//     ├── pol_qc.json             — QC summary (optional; standalone diagnostic
//     │                             file only -- not read by any downstream
//     │                             module or report, browse it directly)
//     ├── README_pol_metrics.txt  — Documentation
//     └── pol_metrics.log         — Processing log
//
// Parameters (params.pol.*, except mapq/dedup which share params.qc.* with
// module 13 so one setting controls MAPQ/dedup filtering everywhere):
//   tss_win            : TSS window ±bp (default: 50)
//   body_offset_min    : Min body start offset (default: 2000 bp)
//   body_offset_frac   : Body start as fraction (default: 0.10)
//   min_body_frac      : Min body window as fraction of gene length for a
//                        valid pausing index; shorter → truncated/NaN (default: 0.10)
//   min_body_len       : Optional absolute body-length floor in bp (default: 0 = off)
//   feature_types      : GTF features to use (default: "gene,transcript")
//   fail_if_no_genes   : Fail on empty gene list (default: false)
//
// ============================================================================


process calculate_polymerase_occupancy_metrics {

  tag        { sid }
  label      'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/08_pol_metrics/${sid}" },
             mode: params.publish_mode,
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sid),
          path(in_bam),
          path(func_bed),
          path(pos3_cpm_bg,   stageAs: 'track_cpm_pos.bedgraph'),
          path(neg3_cpm_bg,   stageAs: 'track_cpm_neg.bedgraph'),
          path(pos3_sicpm_bg, stageAs: 'track_sicpm_pos.bedgraph'),
          path(neg3_sicpm_bg, stageAs: 'track_sicpm_neg.bedgraph'),
          val(cond), val(tp), val(rep)
    // Gene model: the gtf_to_catalog genes.tsv catalog (NOT the raw GTF). Using
    // the same catalog as functional-region calling (module 10) makes the TSS/TES
    // for the pausing index identical to those used for promoter signal.
    path genes_cat

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sid), path('pol_gene_metrics.tsv'), 
          val(cond), val(tp), val(rep), emit: genes
    tuple val(sid), path('pausing_index.tsv'),
          val(cond), val(tp), val(rep), emit: pausing
    tuple val(sid), path('pol_density.tsv'),
          val(cond), val(tp), val(rep), emit: density
    path 'pol_qc.json', optional: true, emit: qc
    path 'README_pol_metrics.txt'
    path 'pol_metrics.log', emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Limit threading for linear algebra libraries
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export BLIS_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a pol_metrics.log)
  exec 2> >(tee -a pol_metrics.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "calculate_polymerase_occupancy_metrics" "Unexpected process failure" "Check pol_metrics.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "POL | START | sample=${sid} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="${sid}"
  CONDITION="${cond}"
  TIMEPOINT="${tp}"
  REPLICATE="${rep}"
  THREADS=${task.cpus}

  # Input files
  IN_BAM="${in_bam}"
  FUNC_BED="${func_bed}"
  GENES_CAT="${genes_cat}"
  CALC_SCRIPT="\$(command -v calculate_pol_metrics.py)"

  # Which mate carries the Pol II 3'-end signal for PE BAMs -- same knob and
  # default as module 06's coverage tracks (params.align.pe_signal_mate), so
  # pausing-index counting stays consistent with the tracks it's compared
  # against. Ignored by the script for single-end BAMs.
  PE_SIGNAL_MATE="${params.align?.pe_signal_mate ?: 'read2'}"

  # Coverage tracks
  POS_CPM="${pos3_cpm_bg}"
  NEG_CPM="${neg3_cpm_bg}"
  POS_SICPM="${pos3_sicpm_bg}"
  NEG_SICPM="${neg3_sicpm_bg}"

  # Parameters
  MAPQ=${params.qc?.mapq ?: 10}
  MULTIMAP_K=${params.align?.multimap_k ?: 0}
  # Uniqueness filter for gene quantification. With bowtie2 -k, MAPQ is set to
  # 255 (unavailable), so uniquely-mapped reads are selected via the NH tag
  # (NH==1, written by add_nh_tags.awk in module 05). In legacy single-best
  # mode there is no NH tag, so fall back to the MAPQ threshold.
  if [[ "\${MULTIMAP_K}" -gt 1 ]]; then
    UNIQUE_FILTER="-d NH:1"
    UNIQUE_DESC="NH==1 (bowtie2 -k mode)"
  else
    UNIQUE_FILTER="-q \${MAPQ}"
    UNIQUE_DESC="MAPQ≥\${MAPQ}"
  fi
  DEDUP_ENABLED=\$([[ "${params.qc?.dedup ?: true}" == "false" ]] && echo 0 || echo 1)
  TSS_WIN=${params.pol?.tss_win ?: 50}
  BODY_OFFSET_MIN=${params.pol?.body_offset_min ?: 2000}
  BODY_OFFSET_FRAC=${params.pol?.body_offset_frac ?: 0.10}
  MIN_BODY_FRAC=${params.pol?.min_body_frac ?: 0.10}
  MIN_BODY_LEN=${params.pol?.min_body_len ?: 0}
  FEATURE_TYPES="${params.pol?.feature_types ?: 'gene,transcript'}"
  FAIL_IF_NO_GENES=\$([[ "${params.pol?.fail_if_no_genes}" == "true" ]] && echo 1 || echo 0)

  echo "POL | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "POL | CONFIG | Condition: \${CONDITION}"
  echo "POL | CONFIG | Timepoint: \${TIMEPOINT}"
  echo "POL | CONFIG | Replicate: \${REPLICATE}"
  echo "POL | CONFIG | Threads: \${THREADS}"
  echo ""
  echo "POL | CONFIG | Input Files:"
  echo "POL | CONFIG |   BAM: \$(basename \${IN_BAM})"
  echo "POL | CONFIG |   Functional regions: \$(basename \${FUNC_BED})"
  echo "POL | CONFIG |   Gene catalog: \$(basename \${GENES_CAT})"
  echo ""
  echo "POL | CONFIG | Coverage Tracks:"
  echo "POL | CONFIG |   CPM: \$(basename \${POS_CPM}), \$(basename \${NEG_CPM})"
  echo "POL | CONFIG |   siCPM: \$(basename \${POS_SICPM}), \$(basename \${NEG_SICPM})"
  echo ""
  echo "POL | CONFIG | BAM Filtering:"
  echo "POL | CONFIG |   MAPQ threshold: \${MAPQ}"
  echo "POL | CONFIG |   Remove duplicates: \$([ \${DEDUP_ENABLED} -eq 1 ] && echo "yes" || echo "no")"
  echo ""
  echo "POL | CONFIG | Gene Metrics:"
  echo "POL | CONFIG |   TSS window: ±\${TSS_WIN} bp"
  echo "POL | CONFIG |   Body offset min: \${BODY_OFFSET_MIN} bp"
  echo "POL | CONFIG |   Body offset fraction: \${BODY_OFFSET_FRAC}"
  echo "POL | CONFIG |   Feature types: \${FEATURE_TYPES}"
  echo "POL | CONFIG |   Fail if no genes: \$([ \${FAIL_IF_NO_GENES} -eq 1 ] && echo "yes" || echo "no")"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "POL | VALIDATE | Checking input files..."

  # Shared resolver (bin/tracktx_error_fragment.sh): micromamba (container) ->
  # /opt/conda (container fallback) -> bare python3 (conda profile/local)
  tracktx_resolve_python

  # Check Python script
  if [[ ! -f "\${CALC_SCRIPT}" ]]; then
    tracktx_error "calculate_polymerase_occupancy_metrics" "Python script not found: \${CALC_SCRIPT}" "Ensure bin/calculate_pol_metrics.py exists"
  fi
  echo "POL | VALIDATE | Python script: \${CALC_SCRIPT}"

  # Check BAM
  if [[ ! -s "\${IN_BAM}" ]]; then
    tracktx_error "calculate_polymerase_occupancy_metrics" "BAM file missing or empty: \${IN_BAM}" "Check upstream alignment module"
  fi
  BAM_SIZE=\$(tracktx_size "\${IN_BAM}")
  echo "POL | VALIDATE | BAM: \${BAM_SIZE} bytes"

  # Check gene catalog
  if [[ ! -s "\${GENES_CAT}" ]]; then
    tracktx_error "calculate_polymerase_occupancy_metrics" "Gene catalog missing or empty: \${GENES_CAT}" "Check download_genome_annotations (genes.tsv) module"
  fi
  GENES_SIZE=\$(tracktx_size "\${GENES_CAT}")
  GENES_LINES=\$(wc -l < "\${GENES_CAT}" 2>/dev/null | tr -d ' ' || echo 0)
  echo "POL | VALIDATE | Gene catalog: \${GENES_SIZE} bytes (\${GENES_LINES} lines)"

  # Check CPM tracks (required)
  for TRACK in "\${POS_CPM}" "\${NEG_CPM}"; do
    if [[ ! -e "\${TRACK}" ]]; then
      tracktx_error "calculate_polymerase_occupancy_metrics" "Required CPM track missing: \${TRACK}" "Check normalize_coverage_tracks module"
    fi
    TRACK_SIZE=\$(tracktx_size "\${TRACK}")
    echo "POL | VALIDATE | \$(basename \${TRACK}): \${TRACK_SIZE} bytes"
  done

  # Check siCPM tracks (optional)
  if [[ -s "\${POS_SICPM}" && -s "\${NEG_SICPM}" ]]; then
    SICPM_SIZE=\$(tracktx_size "\${POS_SICPM}")
    echo "POL | VALIDATE | siCPM tracks available: \${SICPM_SIZE} bytes"
    SICPM_AVAILABLE=1
  else
    echo "POL | VALIDATE | siCPM tracks not available, will use CPM"
    SICPM_AVAILABLE=0
  fi

  # Check functional regions (optional but expected)
  if [[ "\${FUNC_BED}" != "-" && -s "\${FUNC_BED}" ]]; then
    FUNC_SIZE=\$(tracktx_size "\${FUNC_BED}")
    FUNC_COUNT=\$(grep -v '^#' "\${FUNC_BED}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "POL | VALIDATE | Functional regions: \${FUNC_COUNT} regions (\${FUNC_SIZE} bytes)"
  else
    echo "POL | VALIDATE | WARNING: No functional regions provided"
  fi

  # Validate tools
  for TOOL in samtools bedtools awk; do
    if command -v \${TOOL} >/dev/null 2>&1; then
      echo "POL | VALIDATE | \${TOOL}: \$(which \${TOOL})"
    else
      tracktx_error "calculate_polymerase_occupancy_metrics" "Required tool not found: \${TOOL}" "Install \${TOOL} or use -profile docker"
    fi
  done
  if \${PYTHON_CMD} --version >/dev/null 2>&1; then
    echo "POL | VALIDATE | python: \$(\${PYTHON_CMD} --version 2>&1)"
  else
    tracktx_error "calculate_polymerase_occupancy_metrics" "Python not found (tried: \${PYTHON_CMD})" "Use -profile docker"
  fi

  ###########################################################################
  # 3) DENSITY METRICS FROM NORMALIZED TRACKS
  ###########################################################################

  echo "POL | DENSITY | Calculating density metrics from normalized tracks..."

  # Choose between siCPM (preferred) and CPM
  select_track() {
    local sicpm="\$1"
    local cpm="\$2"
    
    if [[ -s "\${sicpm}" ]]; then
      echo "\${sicpm}"
      return 0
    else
      echo "\${cpm}"
      return 0
    fi
  }

  POS_TRACK=\$(select_track "\${POS_SICPM}" "\${POS_CPM}")
  NEG_TRACK=\$(select_track "\${NEG_SICPM}" "\${NEG_CPM}")

  if [[ "\${POS_TRACK}" == "\${POS_SICPM}" ]]; then
    NORM_METHOD="siCPM"
    echo "POL | DENSITY | Using siCPM normalization (spike-in)"
  else
    NORM_METHOD="CPM"
    echo "POL | DENSITY | Using CPM normalization (standard)"
  fi

  # Helper to read bedGraph (handles gzip)
  read_bedgraph() {
    local file="\$1"
    
    if [[ "\${file}" == *.gz ]]; then
      gzip -cd "\${file}"
    else
      cat "\${file}"
    fi
  }

  # Helper to clean and take absolute values
  clean_and_abs() {
    local input="\$1"
    local output="\$2"
    
    if [[ ! -s "\${input}" ]]; then
      : > "\${output}"
      return 0
    fi
    
    echo "POL | DENSITY | Cleaning: \$(basename \${input})"
    
    read_bedgraph "\${input}" | \\
      awk 'BEGIN{OFS="\\t"}
           /^#/ || /^track/ || /^browser/ {next}
           (NF>=4) {
             chr=\$1
             start=\$2+0
             end=\$3+0
             val=\$4+0
             if (end > start) {
               if (val < 0) val = -val
               print chr, start, end, val
             }
           }' | \\
      LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "\${output}"
    
    LINE_COUNT=\$(wc -l < "\${output}" | tr -d ' ')
    echo "POL | DENSITY | Cleaned: \${LINE_COUNT} intervals"
  }

  # Clean and take absolute values of both strands
  clean_and_abs "\${POS_TRACK}" "pos.abs.bedgraph"
  clean_and_abs "\${NEG_TRACK}" "neg.abs.bedgraph"

  # Merge positive and negative strands (|pos| + |neg|)
  echo "POL | DENSITY | Merging strands..."
  
  cat pos.abs.bedgraph neg.abs.bedgraph | \\
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n | \\
    bedtools merge -i - -c 4 -o sum > combined.norm.bedgraph || \\
    : > combined.norm.bedgraph

  COMBINED_LINES=\$(wc -l < combined.norm.bedgraph | tr -d ' ')
  COMBINED_SIZE=\$(tracktx_size combined.norm.bedgraph)
  echo "POL | DENSITY | Combined track: \${COMBINED_LINES} intervals (\${COMBINED_SIZE} bytes)"

  # Map signal to functional regions
  if [[ "\${FUNC_BED}" != "-" && -s "\${FUNC_BED}" ]]; then
    echo "POL | DENSITY | Mapping signal to functional regions..."
    
    # Clean functional regions BED
    awk 'BEGIN{OFS="\\t"}
         !/^track/ && !/^browser/ && !/^#/ && (NF>=3) {
           print \$1, \$2, \$3, (NF>=4 ? \$4 : "."), (NF>=5 ? \$5 : "0"), (NF>=6 ? \$6 : ".")
         }' "\${FUNC_BED}" | \\
      LC_ALL=C sort -k1,1 -k2,2n -k3,3n > functional_regions.sorted.bed
    
    FUNC_SORTED=\$(wc -l < functional_regions.sorted.bed | tr -d ' ')
    echo "POL | DENSITY | Sorted functional regions: \${FUNC_SORTED}"
    
    # Create header
    echo -e "chr\\tstart\\tend\\tname\\tsignal\\tnorm_method" > pol_density.tsv
    
    # Map signal using bedtools
    bedtools map \\
      -a functional_regions.sorted.bed \\
      -b combined.norm.bedgraph \\
      -c 4 \\
      -o sum \\
      -null 0 | \\
      awk -v OFS='\\t' -v METHOD="\${NORM_METHOD}" '{
        print \$1, \$2, \$3, (\$4 != "." ? \$4 : "."), (\$NF != "." ? \$NF : 0), METHOD
      }' >> pol_density.tsv
    
    DENSITY_LINES=\$(tail -n +2 pol_density.tsv | wc -l | tr -d ' ')
    echo "POL | DENSITY | Density table: \${DENSITY_LINES} regions"
  else
    echo "POL | DENSITY | No functional regions, creating header-only file"
    echo -e "chr\\tstart\\tend\\tname\\tsignal\\tnorm_method" > pol_density.tsv
  fi

  ###########################################################################
  # 4) PREPARE BAM FOR GENE METRICS
  ###########################################################################

  echo "POL | BAM | Preparing BAM for gene-level metrics..."

  # Check if BAM is coordinate sorted
  SO_COORD=0
  if samtools view -H "\${IN_BAM}" | \\
     awk '/^@HD/ && /SO:coordinate/ {ok=1} END{exit ok?0:1}'; then
    SO_COORD=1
    echo "POL | BAM | BAM is coordinate-sorted"
  else
    echo "POL | BAM | BAM is not coordinate-sorted, will sort"
  fi

  # Build filtering flags as a SINGLE combined -F value. samtools only honors
  # the LAST -F given if it is passed twice (it does not OR repeated -F
  # values together), so unmapped + duplicate exclusion must be combined into
  # one flag rather than appended as a second -F argument.
  if [[ \${DEDUP_ENABLED} -eq 1 ]]; then
    FILTER_FLAGS="-F 0x404"  # Exclude unmapped (0x4) + duplicates (0x400)
    echo "POL | BAM | Will exclude duplicates"
  else
    FILTER_FLAGS="-F 0x4"  # Exclude unmapped only
    echo "POL | BAM | Will retain duplicates"
  fi

  # Filter BAM
  BAM_START=\$(date +%s)
  
  if [[ \${SO_COORD} -eq 1 ]]; then
    echo "POL | BAM | Filtering BAM (unique=\${UNIQUE_DESC})..."
    samtools view \\
      -@ \${THREADS} \\
      -b \\
      \${UNIQUE_FILTER} \\
      \${FILTER_FLAGS} \\
      "\${IN_BAM}" \\
      -o filtered.bam
  else
    echo "POL | BAM | Filtering and sorting BAM (unique=\${UNIQUE_DESC})..."
    samtools view \\
      -@ \${THREADS} \\
      -b \\
      \${UNIQUE_FILTER} \\
      \${FILTER_FLAGS} \\
      "\${IN_BAM}" | \\
    samtools sort -@ \${THREADS} -o filtered.bam
  fi

  BAM_END=\$(date +%s)
  BAM_TIME=\$((BAM_END - BAM_START))
  
  echo "POL | BAM | Indexing filtered BAM..."
  samtools index -@ \${THREADS} filtered.bam

  FILT_SIZE=\$(tracktx_size filtered.bam)
  FILT_READS=\$(samtools view -c filtered.bam)
  
  echo "POL | BAM | Filtered BAM: \${FILT_SIZE} bytes (\${FILT_READS} reads)"
  echo "POL | BAM | Processing time: \${BAM_TIME}s"

  ###########################################################################
  # 5) CALCULATE GENE METRICS AND PAUSING INDEX
  ###########################################################################

  echo "POL | GENES | Calculating per-gene metrics..."
  echo "POL | GENES | This may take several minutes for large gene sets..."

  GENES_START=\$(date +%s)

  # Build optional flags
  FAIL_FLAG=""
  [[ \${FAIL_IF_NO_GENES} -eq 1 ]] && FAIL_FLAG="--fail-if-empty true"

  set +e
  \${PYTHON_CMD} "\${CALC_SCRIPT}" \\
    --bam filtered.bam \\
    --genes "\${GENES_CAT}" \\
    --pe-signal-mate "\${PE_SIGNAL_MATE}" \\
    --tss-win \${TSS_WIN} \\
    --body-offset-min \${BODY_OFFSET_MIN} \\
    --body-offset-frac \${BODY_OFFSET_FRAC} \\
    --min-body-frac \${MIN_BODY_FRAC} \\
    --min-body-len \${MIN_BODY_LEN} \\
    --feature-types "\${FEATURE_TYPES}" \\
    --out-pausing pausing_index.tsv \\
    --out-genes pol_gene_metrics.tsv \\
    --out-qc pol_qc.json \\
    --threads \${THREADS} \\
    \${FAIL_FLAG}
  
  GENES_RC=\$?
  set -e

  GENES_END=\$(date +%s)
  GENES_TIME=\$((GENES_END - GENES_START))

  echo "POL | GENES | Processing completed in \${GENES_TIME}s"

  # Handle failures
  if [[ \${GENES_RC} -ne 0 ]]; then
    if [[ \${FAIL_IF_NO_GENES} -eq 1 ]]; then
      tracktx_error "calculate_polymerase_occupancy_metrics" "Gene metrics calculation failed with exit code \${GENES_RC} (fail_if_no_genes=true)" "Check pol_metrics.log in work dir" \${GENES_RC}
    else
      echo "POL | WARNING | Creating empty output files"
    fi
  fi

  # Ensure output files exist
  if [[ ! -s pol_gene_metrics.tsv ]]; then
    cat > pol_gene_metrics.tsv <<'GENESEOF'
gene_id	gene_name	chrom	strand	tss_lo	tss_hi	tss_width	body_lo	body_hi	body_len	tss_count	tss_cpm	tss_density_per_bp	body_count	body_cpm	body_density_per_bp	pi_raw	pi_len_norm	is_truncated
GENESEOF
  fi

  if [[ ! -s pausing_index.tsv ]]; then
    cat > pausing_index.tsv <<PAUSINGEOF
gene_id	chrom	strand	tss_count	gene_body_count	pi_raw	pi_len_norm	is_truncated
PAUSINGEOF
  fi

  # Count results
  GENE_COUNT=\$(tail -n +2 pol_gene_metrics.tsv 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  PAUSING_COUNT=\$(tail -n +2 pausing_index.tsv 2>/dev/null | wc -l | tr -d ' ' || echo 0)

  echo "POL | GENES | Gene metrics: \${GENE_COUNT} genes"
  echo "POL | GENES | Pausing indices: \${PAUSING_COUNT} genes"

  ###########################################################################
  # 6) CREATE README
  ###########################################################################

  echo "POL | README | Creating documentation..."

  cat > README_pol_metrics.txt <<DOCEOF
POL-II METRICS — ${sid}
────────────────────────────────────────────────────────────────────────────
  Three views of transcription:
    1. Density: |pos|+|neg| signal (\${NORM_METHOD}, siCPM preferred over CPM
       when available) mapped to functional regions -> pol_density.tsv
       (\$([ -s pol_density.tsv ] && echo "\$(tail -n +2 pol_density.tsv | wc -l) regions" || echo "no data") this run)
    2. Gene metrics: read counts in TSS (+/-\${TSS_WIN}bp) and gene-body windows
       from the BAM (MAPQ>=\${MAPQ}, duplicates \$([ \${DEDUP_ENABLED} -eq 1 ] && echo "removed" || echo "retained"),
       \${FILT_READS} filtered reads) -> pol_gene_metrics.tsv (\${GENE_COUNT} genes).
       Body start = max(\${BODY_OFFSET_MIN}bp, \${BODY_OFFSET_FRAC} x gene_length).
    3. Pausing index: PI = (TSS_density/TSS_width) / (Body_density/Body_width).
       PI>1.5 = strong promoter-proximal pausing; PI<1.0 = productive elongation.
       pi_len_norm is length-normalized (use for cross-gene comparison); pi_raw
       is the direct ratio. is_truncated flags genes whose body window couldn't
       be properly defined. -> pausing_index.tsv (\${PAUSING_COUNT} genes)

  pol_qc.json — standalone diagnostic (not read by any downstream report)

  Feature types: \${FEATURE_TYPES}  |  BAM filter+gene-metric time: \$((BAM_TIME + GENES_TIME))s
DOCEOF

  echo "POL | README | Documentation created"

  ###########################################################################
  # 7) VALIDATION
  ###########################################################################

  echo "POL | VALIDATE | Validating outputs..."

  # Check required outputs
  for FILE in pol_gene_metrics.tsv pausing_index.tsv pol_density.tsv; do
    if [[ ! -e "\${FILE}" ]]; then
      tracktx_error "calculate_polymerase_occupancy_metrics" "Missing output file: \${FILE}" "Check pol_metrics.log in work dir"
    fi
  done

  # Validate file formats
  if [[ -s pol_gene_metrics.tsv ]]; then
    GENES_COLS=\$(head -1 pol_gene_metrics.tsv | awk -F'\\t' '{print NF}')
    if [[ \${GENES_COLS} -ne 19 ]]; then
      echo "POL | WARNING | Gene metrics has \${GENES_COLS} columns, expected 19"
    fi
  fi

  if [[ -s pausing_index.tsv ]]; then
    PAUSING_COLS=\$(head -1 pausing_index.tsv | awk -F'\\t' '{print NF}')
    # calculate_pol_metrics.py emits 8 columns:
    #   gene_id chrom strand tss_count gene_body_count pi_raw pi_len_norm is_truncated
    if [[ \${PAUSING_COLS} -ne 8 ]]; then
      echo "POL | WARNING | Pausing index has \${PAUSING_COLS} columns, expected 8"
    fi
  fi

  if [[ -s pol_density.tsv ]]; then
    DENSITY_COLS=\$(head -1 pol_density.tsv | awk -F'\\t' '{print NF}')
    if [[ \${DENSITY_COLS} -ne 6 ]]; then
      echo "POL | WARNING | Density has \${DENSITY_COLS} columns, expected 6"
    fi
  fi

  echo "POL | VALIDATE | All outputs validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  # Calculate some statistics if possible
  if [[ \${GENE_COUNT} -gt 0 && -s pausing_index.tsv ]]; then
    # Calculate median PI (excluding header and truncated)
    # Use pi_len_norm (col 7) if 8-column format, else pi_raw (col 6); exclude truncated (last col)
    MEDIAN_PI=\$(tail -n +2 pausing_index.tsv | \\
                awk -F'\\t' '\$NF!="1" && \$NF!="True" && ((NF>=8 && \$7!="NA" && \$7+0>0) || (NF==7 && \$6!="NA" && \$6+0>0)) {print (NF>=8?\$7:\$6)}' | \\
                sort -n | \\
                awk '{a[NR]=\$1} END{print (NR%2==1)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}' || echo "NA")
  else
    MEDIAN_PI="NA"
  fi

  echo "────────────────────────────────────────────────────────────────────────"
  echo "POL | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "POL | SUMMARY | Normalization: \${NORM_METHOD}"
  echo "POL | SUMMARY | Filtered reads: \${FILT_READS}"
  echo "POL | SUMMARY | Genes analyzed: \${GENE_COUNT}"
  echo "POL | SUMMARY | Pausing indices: \${PAUSING_COUNT}"
  echo "POL | SUMMARY | Median PI: \${MEDIAN_PI}"
  echo "POL | SUMMARY | Density regions: \$([ -s pol_density.tsv ] && tail -n +2 pol_density.tsv | wc -l || echo 0)"
  echo "POL | SUMMARY | Processing time: \$((BAM_TIME + GENES_TIME))s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "POL | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}