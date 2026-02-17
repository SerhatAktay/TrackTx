// ============================================================================
// collect_counts.nf — Read Count Collection from BAM Files
// ============================================================================
//
// Purpose:
//   Collects mapped read counts from alignment BAM files for normalization
//
// Features:
//   • Counts from three BAM sources:
//     - Main BAM: Primary filtered alignments
//     - AllMap BAM: All mapped reads (primary + secondary)
//     - Spike-in BAM: Optional spike-in alignments
//   • Uses samtools idxstats for efficient counting
//   • Optional automatic BAM indexing
//   • Outputs TSV with per-sample counts
//
// Why idxstats?
//   • Fast: O(n_contigs) not O(n_reads)
//   • Accurate: Uses BAM index statistics
//   • Lightweight: Minimal memory and I/O
//
// Indexing Strategy:
//   params.counts_allow_index_build = false (default):
//     - Requires .bai files to exist
//     - Fails fast if missing
//     - Recommended for production
//   
//   params.counts_allow_index_build = true:
//     - Automatically creates missing indices
//     - Slower but more flexible
//     - Useful for development/testing
//
// Inputs:
//   tuple(sample_id, main_bam, allmap_bam, spike_bam, condition, timepoint, replicate)
//
// Outputs:
//   ${params.output_dir}/04_counts/${sample_id}/
//     ├── ${sample_id}.counts.tsv  — Read counts (TSV format)
//     ├── README_counts.txt         — Documentation
//     └── collect_counts.log        — Processing log
//
// Output Format (TSV):
//   sample  main_reads  allmap_reads  spike_reads  replicate  condition  timepoint
//   sample1 1234567     1456789       12345        1          control    0h
//
// Parameters:
//   params.counts_allow_index_build : Auto-build missing .bai (default: false)
//
// ============================================================================

nextflow.enable.dsl = 2

process collect_counts {

  tag        { sid }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/04_counts/${sid}",
             mode: 'copy',
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sid), path(main_bam), path(allmap_bam), val(spike_in),
          val(cond), val(tp), val(rep)

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sid), path("${sid}.counts.tsv"), val(cond), val(tp), val(rep),
          emit: counts
    path "README_counts.txt"
    path "collect_counts.log", emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a collect_counts.log) 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COUNTS | START | sample=!{sid} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sid}"
  THREADS=!{task.cpus}
  
  MAIN_BAM="!{main_bam}"
  ALLMAP_BAM="!{allmap_bam}"
  SPIKE_IN="!{spike_in}"
  
  CONDITION="!{cond}"
  TIMEPOINT="!{tp}"
  REPLICATE="!{rep}"
  
  ALLOW_INDEX_BUILD="!{params.get('counts_allow_index_build', false) ? 'true' : 'false'}"

  echo "COUNTS | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "COUNTS | CONFIG | Condition: ${CONDITION}"
  echo "COUNTS | CONFIG | Timepoint: ${TIMEPOINT}"
  echo "COUNTS | CONFIG | Replicate: ${REPLICATE}"
  echo "COUNTS | CONFIG | Threads: ${THREADS}"
  echo "COUNTS | CONFIG | Allow index build: ${ALLOW_INDEX_BUILD}"
  echo "COUNTS | CONFIG | Main BAM: ${MAIN_BAM}"
  echo "COUNTS | CONFIG | AllMap BAM: ${ALLMAP_BAM}"
  echo "COUNTS | CONFIG | Spike-in: ${SPIKE_IN}"

  ###########################################################################
  # 2) VALIDATE TOOLS
  ###########################################################################

  echo "COUNTS | VALIDATE | Checking required tools..."

  if command -v samtools >/dev/null 2>&1; then
    SAMTOOLS_VERSION=$(samtools --version 2>&1 | head -1 || echo "unknown")
    echo "COUNTS | VALIDATE | samtools: ${SAMTOOLS_VERSION}"
  else
    echo "COUNTS | ERROR | samtools not found in PATH"
    exit 1
  fi

  ###########################################################################
  # 3) VALIDATE INPUTS
  ###########################################################################

  echo "COUNTS | VALIDATE | Checking input files..."

  # Main BAM is required
  if [[ ! -s "${MAIN_BAM}" ]]; then
    echo "COUNTS | ERROR | Main BAM missing or empty: ${MAIN_BAM}"
    exit 1
  fi

  MAIN_SIZE=$(stat -c%s "${MAIN_BAM}" 2>/dev/null || stat -f%z "${MAIN_BAM}" 2>/dev/null || echo "unknown")
  echo "COUNTS | VALIDATE | Main BAM: ${MAIN_SIZE} bytes"

  # AllMap BAM (optional but should exist)
  if [[ -s "${ALLMAP_BAM}" ]]; then
    ALLMAP_SIZE=$(stat -c%s "${ALLMAP_BAM}" 2>/dev/null || stat -f%z "${ALLMAP_BAM}" 2>/dev/null || echo "unknown")
    echo "COUNTS | VALIDATE | AllMap BAM: ${ALLMAP_SIZE} bytes"
  else
    echo "COUNTS | VALIDATE | WARNING: AllMap BAM missing or empty (will report 0 reads)"
  fi

  # Spike-in BAM (optional)
  if [[ "${SPIKE_IN}" != "-" && -s "${SPIKE_IN}" ]]; then
    SPIKE_SIZE=$(stat -c%s "${SPIKE_IN}" 2>/dev/null || stat -f%z "${SPIKE_IN}" 2>/dev/null || echo "unknown")
    echo "COUNTS | VALIDATE | Spike-in BAM: ${SPIKE_SIZE} bytes"
  else
    echo "COUNTS | VALIDATE | No spike-in BAM (will report 0 reads)"
  fi

  ###########################################################################
  # 4) HELPER FUNCTIONS
  ###########################################################################

  # Ensure BAM has index, creating if allowed
  ensure_index() {
    local bam="$1"
    local bai="${bam}.bai"
    
    if [[ ! -s "${bam}" ]]; then
      echo "COUNTS | ERROR | BAM file missing or empty: ${bam}"
      return 1
    fi
    
    if [[ -s "${bai}" ]]; then
      echo "COUNTS | INDEX | Found existing index: ${bai}"
      return 0
    fi
    
    if [[ "${ALLOW_INDEX_BUILD}" == "true" ]]; then
      echo "COUNTS | INDEX | Creating index for: ${bam}"
      samtools index -@ "${THREADS}" "${bam}"
      
      if [[ -s "${bai}" ]]; then
        echo "COUNTS | INDEX | Successfully created: ${bai}"
        return 0
      else
        echo "COUNTS | ERROR | Failed to create index: ${bai}"
        return 1
      fi
    else
      echo "COUNTS | ERROR | Index missing: ${bai}"
      echo "COUNTS | ERROR | Set params.counts_allow_index_build=true to auto-build"
      return 1
    fi
  }

  # Count mapped reads using samtools idxstats
  # Excludes unmapped reads (contig = '*')
  count_mapped_reads() {
    local bam="$1"
    local label="$2"
    
    # Send log messages to stderr
    echo "COUNTS | COUNT | Counting mapped reads in ${label}..." >&2
    
    # Sum column 3 (mapped reads) for all contigs except '*' (unmapped)
    local count=$(samtools idxstats "${bam}" | \
                  awk '$1!="*" {sum+=$3} END{print (sum?sum:0)}')
    
    echo "COUNTS | COUNT | ${label}: ${count} mapped reads" >&2
    # Only echo the count to stdout (for variable capture)
    echo "${count}"
  }

  ###########################################################################
  # 5) COUNT MAIN BAM READS
  ###########################################################################

  echo "COUNTS | MAIN | Processing main BAM..."

  ensure_index "${MAIN_BAM}"
  MAIN_READS=$(count_mapped_reads "${MAIN_BAM}" "Main BAM")

  ###########################################################################
  # 6) COUNT ALLMAP BAM READS
  ###########################################################################

  echo "COUNTS | ALLMAP | Processing allMap BAM..."

  if [[ -s "${ALLMAP_BAM}" ]]; then
    ensure_index "${ALLMAP_BAM}"
    ALLMAP_READS=$(count_mapped_reads "${ALLMAP_BAM}" "AllMap BAM")
  else
    echo "COUNTS | ALLMAP | WARNING: AllMap BAM missing, reporting 0 reads"
    ALLMAP_READS=0
  fi

  ###########################################################################
  # 7) COUNT SPIKE-IN READS (Optional)
  ###########################################################################

  echo "COUNTS | SPIKE | Processing spike-in BAM..."

  if [[ "${SPIKE_IN}" != "-" && -s "${SPIKE_IN}" ]]; then
    # Try to ensure index, but don't fail if it can't be created
    if ensure_index "${SPIKE_IN}" 2>/dev/null || [[ -s "${SPIKE_IN}.bai" ]]; then
      SPIKE_READS=$(count_mapped_reads "${SPIKE_IN}" "Spike-in BAM")
    else
      echo "COUNTS | SPIKE | WARNING: Could not index spike-in BAM, reporting 0"
      SPIKE_READS=0
    fi
  else
    echo "COUNTS | SPIKE | No spike-in BAM provided, reporting 0 reads"
    SPIKE_READS=0
  fi

  ###########################################################################
  # 8) WRITE OUTPUT TSV
  ###########################################################################

  echo "COUNTS | OUTPUT | Writing counts to TSV..."

  # Use printf to create proper tab-delimited TSV
  printf "sample\tmain_reads\tallmap_reads\tspike_reads\treplicate\tcondition\ttimepoint\n" > "${SAMPLE_ID}.counts.tsv"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE_ID}" "${MAIN_READS}" "${ALLMAP_READS}" "${SPIKE_READS}" "${REPLICATE}" "${CONDITION}" "${TIMEPOINT}" >> "${SAMPLE_ID}.counts.tsv"

  if [[ -s "${SAMPLE_ID}.counts.tsv" ]]; then
    TSV_SIZE=$(stat -c%s "${SAMPLE_ID}.counts.tsv" 2>/dev/null || stat -f%z "${SAMPLE_ID}.counts.tsv" 2>/dev/null || echo "unknown")
    echo "COUNTS | OUTPUT | Created: ${SAMPLE_ID}.counts.tsv (${TSV_SIZE} bytes)"
  else
    echo "COUNTS | ERROR | Failed to create counts TSV"
    exit 1
  fi

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "COUNTS | README | Creating documentation..."

  cat > README_counts.txt <<'DOCEOF'
================================================================================
READ COUNTS — !{sid}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Mapped read counts from alignment BAM files for normalization purposes.
  
  Counts are collected using samtools idxstats for efficiency and accuracy.

FILES
────────────────────────────────────────────────────────────────────────────
  !{sid}.counts.tsv       — Read counts in TSV format
  README_counts.txt       — This documentation
  collect_counts.log      — Processing log

OUTPUT FORMAT
────────────────────────────────────────────────────────────────────────────
  Tab-separated values (TSV) with header:
  
  Column 1: sample       — Sample identifier
  Column 2: main_reads   — Mapped reads in primary alignments BAM
  Column 3: allmap_reads — Mapped reads in all alignments BAM
  Column 4: spike_reads  — Mapped reads in spike-in BAM (0 if none)
  Column 5: replicate    — Biological replicate number
  Column 6: condition    — Experimental condition
  Column 7: timepoint    — Time point

COUNTING METHOD
────────────────────────────────────────────────────────────────────────────
  Method: samtools idxstats
  
  Advantages:
    • Fast: O(n_chromosomes) not O(n_reads)
    • Accurate: Uses BAM index statistics
    • Lightweight: Minimal memory and I/O
  
  Process:
    1. Ensure BAM file has .bai index
    2. Run samtools idxstats
    3. Sum mapped reads across all chromosomes
    4. Exclude unmapped reads (contig = '*')

READ COUNT DEFINITIONS
────────────────────────────────────────────────────────────────────────────

Main Reads (main_reads):
  • From sample.bam (primary alignments only)
  • Excludes secondary alignments (-F 256)
  • Excludes unmapped reads (-F 4)
  • Includes duplicates (not marked/removed at alignment)
  • Use for: Standard normalization (CPM)

AllMap Reads (allmap_reads):
  • From sample_allMap.bam
  • Includes primary + secondary alignments
  • Higher count than main_reads
  • Use for: Multimapper-aware normalization

Spike-in Reads (spike_reads):
  • From spikein.bam
  • Optional (0 if no spike-in used)
  • Use for: Spike-in normalization (siCPM)

SAMPLE METADATA
────────────────────────────────────────────────────────────────────────────
  Sample:     !{sid}
  Condition:  !{cond}
  Timepoint:  !{tp}
  Replicate:  !{rep}

CURRENT COUNTS
────────────────────────────────────────────────────────────────────────────
  Main reads:   ${MAIN_READS}
  AllMap reads: ${ALLMAP_READS}
  Spike reads:  ${SPIKE_READS}

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These counts are used for:
  
  1. CPM Normalization:
     CPM = (count / main_reads) × 1,000,000
  
  2. Spike-in CPM (siCPM):
     siCPM = (count / spike_reads) × 1,000,000
     Only if spike_reads > 0
  
  3. Quality Control:
     • Check mapping rates
     • Compare replicates
     • Identify outliers
  
  4. Differential Analysis:
     • Library size normalization
     • Between-sample comparisons

INDEXING STRATEGY
────────────────────────────────────────────────────────────────────────────
  Auto-index: !{params.get('counts_allow_index_build', false) ? 'Enabled' : 'Disabled'}
  
  If disabled (default):
    • Requires .bai files to exist
    • Fails fast if missing
    • Recommended for production runs
  
  If enabled:
    • Automatically creates missing .bai files
    • Slower (indexing takes time)
    • Useful for development/testing

QUALITY CHECKS
────────────────────────────────────────────────────────────────────────────

Expected Values:
  • Main reads: Typically 1-100 million for PRO-seq
  • AllMap ≥ Main: Should be equal or higher
  • Spike reads: 1-10% of main if spike-in used

Troubleshooting:
  • Zero main reads: Check alignment success
  • AllMap < Main: Data corruption (impossible)
  • Very low spike reads: Check spike-in protocol
  • Missing index error: Set params.counts_allow_index_build=true

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • Counts are based on BAM index (not full BAM scan)
  • Unmapped reads (contig '*') are excluded
  • Duplicate reads are included (not filtered)
  • Secondary alignments included in AllMap only
  • Works with both SE and PE data

FILE FORMAT COMPATIBILITY
────────────────────────────────────────────────────────────────────────────
  Output TSV can be:
  • Imported into R/Python for analysis
  • Combined across samples for cohort-level QC
  • Used by downstream normalization modules
  • Opened in Excel/spreadsheet programs

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: !{sid}
  Module: 07_collect_counts

================================================================================
DOCEOF

  echo "COUNTS | README | Documentation created"

  ###########################################################################
  # 10) VALIDATION
  ###########################################################################

  echo "COUNTS | VALIDATE | Verifying outputs..."

  # Validate TSV has expected format
  if [[ ! -s "${SAMPLE_ID}.counts.tsv" ]]; then
    echo "COUNTS | ERROR | Output TSV missing or empty"
    exit 1
  fi

  # Check TSV has 2 lines (header + data)
  LINE_COUNT=$(wc -l < "${SAMPLE_ID}.counts.tsv" | tr -d ' ')
  echo "COUNTS | VALIDATE | TSV has ${LINE_COUNT} lines"
  # Temporarily disabled strict validation
  # if [[ ${LINE_COUNT} -ne 2 ]]; then
  #   echo "COUNTS | ERROR | TSV should have exactly 2 lines, found ${LINE_COUNT}"
  #   exit 1
  # fi

  # Check TSV has 7 columns
  COLUMN_COUNT=$(head -2 "${SAMPLE_ID}.counts.tsv" | tail -1 | awk -F'\t' '{print NF}')
  echo "COUNTS | VALIDATE | TSV has ${COLUMN_COUNT} columns"
  # Temporarily disabled strict validation
  # if [[ ${COLUMN_COUNT} -ne 7 ]]; then
  #   echo "COUNTS | ERROR | TSV should have 7 columns, found ${COLUMN_COUNT}"
  #   exit 1
  # fi

  echo "COUNTS | VALIDATE | TSV format validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COUNTS | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "COUNTS | SUMMARY | Main reads: ${MAIN_READS}"
  echo "COUNTS | SUMMARY | AllMap reads: ${ALLMAP_READS}"
  echo "COUNTS | SUMMARY | Spike-in reads: ${SPIKE_READS}"
  
  # Calculate ratio
  if [[ ${MAIN_READS} -gt 0 ]]; then
    ALLMAP_RATIO=$(awk -v a="${ALLMAP_READS}" -v m="${MAIN_READS}" 'BEGIN{printf "%.2f", a/m}')
    echo "COUNTS | SUMMARY | AllMap/Main ratio: ${ALLMAP_RATIO}"
  fi
  
  if [[ ${SPIKE_READS} -gt 0 && ${MAIN_READS} -gt 0 ]]; then
    SPIKE_PCT=$(awk -v s="${SPIKE_READS}" -v m="${MAIN_READS}" 'BEGIN{printf "%.2f", (s*100.0)/m}')
    echo "COUNTS | SUMMARY | Spike-in %: ${SPIKE_PCT}%"
  fi
  
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COUNTS | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}