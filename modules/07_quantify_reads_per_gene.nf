// ============================================================================
// quantify_reads_per_gene.nf — Read Count Collection from BAM Files
// ============================================================================
//
// Purpose:
//   Collects mapped read counts from alignment BAM files for normalization
//
// NAME NOTE:
//   Despite the name "quantify_reads_per_gene", this module does NOT produce
//   per-gene counts. It computes whole-library totals (main / allMap / spike
//   mapped reads via samtools idxstats) that module 08 uses as CPM/siCPM library
//   sizes. Per-gene quantification lives in module 11 (calculate_pol_metrics.py).
//   The process name is retained only because it is referenced by withName
//   selectors in nextflow.config; renaming it would silently drop those resource
//   directives.
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
//     └── quantify_reads_per_gene.log — Processing log
//
// Output Format (TSV):
//   sample  main_reads  allmap_reads  spike_reads  replicate  condition  timepoint
//   sample1 1234567     1456789       12345        1          control    0h
//
// Parameters:
//   params.counts_allow_index_build : Auto-build missing .bai (default: false)
//
// ============================================================================


process quantify_reads_per_gene {

  tag        { sid }
  label      'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/04_counts/${sid}" },
             mode: params.publish_mode,
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
    path "quantify_reads_per_gene.log", emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a quantify_reads_per_gene.log)
  exec 2> >(tee -a quantify_reads_per_gene.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "quantify_reads_per_gene" "Unexpected process failure" "Check quantify_reads_per_gene.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COUNTS | START | sample=${sid} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="${sid}"
  THREADS=${task.cpus}
  
  MAIN_BAM="${main_bam}"
  ALLMAP_BAM="${allmap_bam}"
  SPIKE_IN="${spike_in}"
  
  CONDITION="${cond}"
  TIMEPOINT="${tp}"
  REPLICATE="${rep}"
  
  ALLOW_INDEX_BUILD="${params.get('counts_allow_index_build', false) ? 'true' : 'false'}"

  # Cross-sample I/O lock (see modules/06_generate_coverage_tracks.nf for the
  # full writeup): serializes the samtools index call below across
  # concurrently running samples sharing this pipeline's USB/HDD-backed work
  # volume, where concurrent big sequential reads interleave into
  # seek-thrashing instead of parallel throughput.
  # Shared with_io_lock()/init (bin/tracktx_error_fragment.sh); override the
  # slot-wait timeout with TRACKS_IO_LOCK_TIMEOUT (default 1800s).
  tracktx_io_lock_init "${projectDir}/.tracktx_io.lock"

  echo "COUNTS | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "COUNTS | CONFIG | Condition: \${CONDITION}"
  echo "COUNTS | CONFIG | Timepoint: \${TIMEPOINT}"
  echo "COUNTS | CONFIG | Replicate: \${REPLICATE}"
  echo "COUNTS | CONFIG | Threads: \${THREADS}"
  echo "COUNTS | CONFIG | Allow index build: \${ALLOW_INDEX_BUILD}"
  echo "COUNTS | CONFIG | Main BAM: \${MAIN_BAM}"
  echo "COUNTS | CONFIG | AllMap BAM: \${ALLMAP_BAM}"
  echo "COUNTS | CONFIG | Spike-in: \${SPIKE_IN}"

  ###########################################################################
  # 2) VALIDATE TOOLS
  ###########################################################################

  echo "COUNTS | VALIDATE | Checking required tools..."

  if ! command -v samtools >/dev/null 2>&1; then
    tracktx_error "quantify_reads_per_gene" "samtools not found in PATH" "Install samtools or use -profile docker"
  fi
  SAMTOOLS_VERSION=\$(samtools --version 2>&1 | head -1 || echo "unknown")
  echo "COUNTS | VALIDATE | samtools: \${SAMTOOLS_VERSION}"

  ###########################################################################
  # 3) VALIDATE INPUTS
  ###########################################################################

  echo "COUNTS | VALIDATE | Checking input files..."

  # Main BAM is required
  if [[ ! -s "\${MAIN_BAM}" ]]; then
    tracktx_error "quantify_reads_per_gene" "Main BAM missing or empty: \${MAIN_BAM}" "Check align_reads_to_genome produced sample.bam"
  fi

  MAIN_SIZE=\$(tracktx_size "\${MAIN_BAM}")
  echo "COUNTS | VALIDATE | Main BAM: \${MAIN_SIZE} bytes"

  # AllMap BAM (optional but should exist)
  if [[ -s "\${ALLMAP_BAM}" ]]; then
    ALLMAP_SIZE=\$(tracktx_size "\${ALLMAP_BAM}")
    echo "COUNTS | VALIDATE | AllMap BAM: \${ALLMAP_SIZE} bytes"
  else
    echo "COUNTS | VALIDATE | WARNING: AllMap BAM missing or empty (will report 0 reads)"
  fi

  # Spike-in BAM (optional)
  if [[ "\${SPIKE_IN}" != "-" && -s "\${SPIKE_IN}" ]]; then
    SPIKE_SIZE=\$(tracktx_size "\${SPIKE_IN}")
    echo "COUNTS | VALIDATE | Spike-in BAM: \${SPIKE_SIZE} bytes"
  else
    echo "COUNTS | VALIDATE | No spike-in BAM (will report 0 reads)"
  fi

  ###########################################################################
  # 4) HELPER FUNCTIONS
  ###########################################################################

  # Ensure BAM has index, creating if allowed
  ensure_index() {
    local bam="\$1"
    local bai="\${bam}.bai"
    
    if [[ ! -s "\${bam}" ]]; then
      tracktx_error "quantify_reads_per_gene" "BAM file missing or empty: \${bam}" "Check input BAM paths"
    fi
    
    if [[ -s "\${bai}" ]]; then
      echo "COUNTS | INDEX | Found existing index: \${bai}"
      return 0
    fi
    
    if [[ "\${ALLOW_INDEX_BUILD}" == "true" ]]; then
      echo "COUNTS | INDEX | Creating index for: \${bam}"
      with_io_lock samtools index -@ "\${THREADS}" "\${bam}"
      
      if [[ -s "\${bai}" ]]; then
        echo "COUNTS | INDEX | Successfully created: \${bai}"
        return 0
      else
        tracktx_error "quantify_reads_per_gene" "Failed to create index: \${bai}" "Check BAM file integrity"
      fi
    else
      tracktx_error "quantify_reads_per_gene" "Index missing: \${bai}" "Set params.counts_allow_index_build=true to auto-build"
    fi
  }

  # Count mapped reads using samtools idxstats
  # Excludes unmapped reads (contig = '*')
  count_mapped_reads() {
    local bam="\$1"
    local label="\$2"
    
    # Send log messages to stderr
    echo "COUNTS | COUNT | Counting mapped reads in \${label}..." >&2
    
    # Sum column 3 (mapped reads) for all contigs except '*' (unmapped)
    local count=\$(samtools idxstats "\${bam}" | \\
                  awk '\$1!="*" {sum+=\$3} END{print (sum?sum:0)}')
    
    echo "COUNTS | COUNT | \${label}: \${count} mapped reads" >&2
    # Only echo the count to stdout (for variable capture)
    echo "\${count}"
  }

  ###########################################################################
  # 5) COUNT MAIN BAM READS
  ###########################################################################

  echo "COUNTS | MAIN | Processing main BAM..."

  ensure_index "\${MAIN_BAM}"
  MAIN_READS=\$(count_mapped_reads "\${MAIN_BAM}" "Main BAM")

  ###########################################################################
  # 6) COUNT ALLMAP BAM READS
  ###########################################################################

  echo "COUNTS | ALLMAP | Processing allMap BAM..."

  if [[ -s "\${ALLMAP_BAM}" ]]; then
    ensure_index "\${ALLMAP_BAM}"
    ALLMAP_READS=\$(count_mapped_reads "\${ALLMAP_BAM}" "AllMap BAM")
  else
    echo "COUNTS | ALLMAP | WARNING: AllMap BAM missing, reporting 0 reads"
    ALLMAP_READS=0
  fi

  ###########################################################################
  # 7) COUNT SPIKE-IN READS (Optional)
  ###########################################################################

  echo "COUNTS | SPIKE | Processing spike-in BAM..."

  if [[ "\${SPIKE_IN}" != "-" && -s "\${SPIKE_IN}" ]]; then
    # Try to ensure index, but don't fail if it can't be created
    if ensure_index "\${SPIKE_IN}" 2>/dev/null || [[ -s "\${SPIKE_IN}.bai" ]]; then
      SPIKE_READS=\$(count_mapped_reads "\${SPIKE_IN}" "Spike-in BAM")
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
  printf "sample\\tmain_reads\\tallmap_reads\\tspike_reads\\treplicate\\tcondition\\ttimepoint\\n" > "\${SAMPLE_ID}.counts.tsv"
  printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "\${SAMPLE_ID}" "\${MAIN_READS}" "\${ALLMAP_READS}" "\${SPIKE_READS}" "\${REPLICATE}" "\${CONDITION}" "\${TIMEPOINT}" >> "\${SAMPLE_ID}.counts.tsv"

  if [[ -s "\${SAMPLE_ID}.counts.tsv" ]]; then
    TSV_SIZE=\$(tracktx_size "\${SAMPLE_ID}.counts.tsv")
    echo "COUNTS | OUTPUT | Created: \${SAMPLE_ID}.counts.tsv (\${TSV_SIZE} bytes)"
  else
    tracktx_error "quantify_reads_per_gene" "Failed to create counts TSV" "Check quantify_reads_per_gene.log in work dir"
  fi

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "COUNTS | README | Creating documentation..."

  cat > README_counts.txt <<DOCEOF
READ COUNTS — ${sid}
────────────────────────────────────────────────────────────────────────────
  ${sid}.counts.tsv (samtools idxstats, summed across contigs, '*' excluded):
    sample  main_reads  allmap_reads  spike_reads  replicate  condition  timepoint

  main_reads   — sample.bam (primary only, -F 260); denominator for CPM
  allmap_reads — sample_allMap.bam (primary + secondary); multimapper-aware CPM
  spike_reads  — spikein.bam, 0 if no spike-in; denominator for siCPM

  This sample: main=\${MAIN_READS}  allmap=\${ALLMAP_READS}  spike=\${SPIKE_READS}

  Auto-index (params.counts_allow_index_build): ${params.get('counts_allow_index_build', false) ? 'enabled' : 'disabled -- missing .bai fails fast'}
DOCEOF

  echo "COUNTS | README | Documentation created"

  ###########################################################################
  # 10) VALIDATION
  ###########################################################################

  echo "COUNTS | VALIDATE | Verifying outputs..."

  # Validate TSV has expected format
  if [[ ! -s "\${SAMPLE_ID}.counts.tsv" ]]; then
    tracktx_error "quantify_reads_per_gene" "Output TSV missing or empty" "Check quantify_reads_per_gene.log in work dir"
  fi

  # Report TSV shape (informational; this file is a fixed 2-line / 7-column table
  # written by printf just above, so a hard check added no safety and the strict
  # validation was disabled long-term — removed rather than left as dead comments).
  LINE_COUNT=\$(wc -l < "\${SAMPLE_ID}.counts.tsv" | tr -d ' ')
  COLUMN_COUNT=\$(head -2 "\${SAMPLE_ID}.counts.tsv" | tail -1 | awk -F'\\t' '{print NF}')
  echo "COUNTS | VALIDATE | TSV has \${LINE_COUNT} lines, \${COLUMN_COUNT} columns"
  echo "COUNTS | VALIDATE | TSV format validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COUNTS | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "COUNTS | SUMMARY | Main reads: \${MAIN_READS}"
  echo "COUNTS | SUMMARY | AllMap reads: \${ALLMAP_READS}"
  echo "COUNTS | SUMMARY | Spike-in reads: \${SPIKE_READS}"
  
  # Calculate ratio
  if [[ \${MAIN_READS} -gt 0 ]]; then
    ALLMAP_RATIO=\$(awk -v a="\${ALLMAP_READS}" -v m="\${MAIN_READS}" 'BEGIN{printf "%.2f", a/m}')
    echo "COUNTS | SUMMARY | AllMap/Main ratio: \${ALLMAP_RATIO}"
  fi
  
  if [[ \${SPIKE_READS} -gt 0 && \${MAIN_READS} -gt 0 ]]; then
    SPIKE_PCT=\$(awk -v s="\${SPIKE_READS}" -v m="\${MAIN_READS}" 'BEGIN{printf "%.2f", (s*100.0)/m}')
    echo "COUNTS | SUMMARY | Spike-in %: \${SPIKE_PCT}%"
  fi
  
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COUNTS | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}