// ============================================================================
// qc_pol_tracktx.nf — Per-Sample Quality Control for PRO-seq Data
// ============================================================================
//
// Purpose:
//   Comprehensive quality control analysis for aligned PRO-seq BAM files
//
// Features:
//   • Alignment statistics (mapping rate, duplicates, MAPQ filtering)
//   • Strand bias analysis (essential for PRO-seq validation)
//   • Fragment length distribution (paired-end only)
//   • Coverage statistics (mean depth across genome)
//   • UMI deduplication statistics (if applicable)
//   • JSON summary for downstream aggregation
//
// Quality Metrics Calculated:
//   1. Mapping Metrics:
//      - Total reads
//      - Mapped reads and rate
//      - Duplicate rate
//      - MAPQ pass rate
//
//   2. Strand Metrics:
//      - Positive strand reads
//      - Negative strand reads
//      - Strand balance (should be ~50/50 for good PRO-seq)
//
//   3. Coverage Metrics:
//      - Mean coverage depth
//      - Per-chromosome coverage
//
//   4. Fragment Metrics (PE only):
//      - Insert size distribution
//      - Expected peak at 150-250 bp
//
//   5. UMI Metrics (if enabled):
//      - Input reads
//      - Deduplicated reads
//      - Duplication rate
//
// Filtering Applied:
//   • Primary alignments only (-F 0x900)
//   • Mapped reads only (-F 0x4)
//   • MAPQ ≥ threshold (default: 10)
//   • Optional duplicate removal (-F 0x400)
//
// Inputs:
//   tuple(sample_id, aligned_bam, dedup_stats,
//         condition, timepoint, replicate)
//
// Outputs:
//   ${params.output_dir}/10_qc/${sample_id}/
//     ├── qc_strand_bias.tsv       — Strand read counts
//     ├── qc_fragment_length.tsv   — Insert size distribution (PE)
//     ├── qc_coverage.tsv          — Coverage statistics
//     ├── qc_pol.json             — JSON summary (all metrics)
//     ├── README_qc.txt            — Documentation
//     └── qc.log                   — Processing log
//
// Parameters (params.qc.*):
//   mapq  : Minimum MAPQ threshold (default: 10)
//   dedup : Remove duplicates (default: true)
//
// Expected Values (PRO-seq):
//   • Mapping rate: >70%
//   • Strand balance: 45-55% (not perfectly balanced due to biology)
//   • Duplicate rate: <30% (varies with depth)
//   • Mean depth: >10× for good coverage
//
// ============================================================================

nextflow.enable.dsl = 2

process qc_pol_tracktx {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/10_qc/${sample_id}",
             mode: 'copy',
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(aligned_bam),
          path(dedup_stats),
          val(condition), val(timepoint), val(replicate)

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path('qc_strand_bias.tsv'),
          path('qc_fragment_length.tsv'),
          path('qc_coverage.tsv'),
          val(condition), val(timepoint), val(replicate), emit: tables
    path 'qc_pol.json',                                     emit: json
    tuple val(sample_id), path('qc_pol.json'),
          val(condition), val(timepoint), val(replicate),    emit: json_meta
    path 'README_qc.txt',                                    emit: readme
    path 'qc.log',                                           emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  def mapq_thr = params.qc?.mapq ?: 10
  def dedup_flag = (params.qc?.dedup == null || params.qc.dedup) ? '-F 0x400' : ''
  def dedup_enabled = (params.qc?.dedup == null || params.qc.dedup)
  
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a qc.log) 2>&1

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "QC | START | sample=${sample_id} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="${sample_id}"
  CONDITION="${condition}"
  TIMEPOINT="${timepoint}"
  REPLICATE="${replicate}"
  
  BAM_FILE="${aligned_bam}"
  DEDUP_STATS="${dedup_stats}"
  
  MAPQ_THRESHOLD=${mapq_thr}
  DEDUP_ENABLED=${dedup_enabled}
  DEDUP_FLAG="${dedup_flag}"

  echo "QC | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "QC | CONFIG | Condition: \${CONDITION}"
  echo "QC | CONFIG | Timepoint: \${TIMEPOINT}"
  echo "QC | CONFIG | Replicate: \${REPLICATE}"
  echo ""
  echo "QC | CONFIG | BAM file: \$(basename \${BAM_FILE})"
  echo "QC | CONFIG | Dedup stats: \$(basename \${DEDUP_STATS})"
  echo ""
  echo "QC | CONFIG | MAPQ threshold: \${MAPQ_THRESHOLD}"
  echo "QC | CONFIG | Remove duplicates: \${DEDUP_ENABLED}"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "QC | VALIDATE | Checking input files..."

  VALIDATION_OK=1

  # Check BAM
  if [[ ! -s "\${BAM_FILE}" ]]; then
    echo "QC | ERROR | BAM file missing or empty: \${BAM_FILE}"
    VALIDATION_OK=0
  else
    BAM_SIZE=\$(stat -c%s "\${BAM_FILE}" 2>/dev/null || stat -f%z "\${BAM_FILE}" 2>/dev/null || echo "unknown")
    echo "QC | VALIDATE | BAM: \${BAM_SIZE} bytes"
  fi

  # Check if BAM is indexed (create if needed)
  if [[ ! -e "\${BAM_FILE}.bai" ]]; then
    echo "QC | VALIDATE | BAM index not found, creating..."
    samtools index "\${BAM_FILE}"
    echo "QC | VALIDATE | BAM index created"
  else
    echo "QC | VALIDATE | BAM index: present"
  fi

  # Check dedup stats (optional)
  if [[ -s "\${DEDUP_STATS}" ]]; then
    DEDUP_SIZE=\$(stat -c%s "\${DEDUP_STATS}" 2>/dev/null || stat -f%z "\${DEDUP_STATS}" 2>/dev/null || echo "unknown")
    echo "QC | VALIDATE | Dedup stats: \${DEDUP_SIZE} bytes"
  else
    echo "QC | VALIDATE | Dedup stats: not available"
  fi

  # Validate tools
  for TOOL in samtools awk; do
    if command -v \${TOOL} >/dev/null 2>&1; then
      echo "QC | VALIDATE | \${TOOL}: \$(which \${TOOL})"
    else
      echo "QC | ERROR | Required tool not found: \${TOOL}"
      VALIDATION_OK=0
    fi
  done

  if [[ \${VALIDATION_OK} -eq 0 ]]; then
    echo "QC | ERROR | Validation failed"
    exit 1
  fi

  ###########################################################################
  # 3) DETECT SEQUENCING MODE
  ###########################################################################

  echo "QC | DETECT | Determining sequencing mode..."

  # Check for paired-end by looking for bowtie2 in BAM header
  IS_PAIRED=0
  if samtools view -H "\${BAM_FILE}" | grep -q '@PG.*bowtie2'; then
    # Additional check: look for properly paired reads
    PAIRED_COUNT=\$(samtools view -c -f 0x1 "\${BAM_FILE}" | head -1)
    if [[ \${PAIRED_COUNT} -gt 0 ]]; then
      IS_PAIRED=1
      echo "QC | DETECT | Sequencing mode: Paired-end"
    else
      echo "QC | DETECT | Sequencing mode: Single-end"
    fi
  else
    echo "QC | DETECT | Sequencing mode: Single-end"
  fi

  ###########################################################################
  # 4) CALCULATE BASIC ALIGNMENT STATISTICS
  ###########################################################################

  echo "QC | STATS | Calculating alignment statistics..."

  STATS_START=\$(date +%s)

  # Total reads (primary alignments, not secondary/supplementary)
  echo "QC | STATS | Counting total reads..."
  TOTAL_READS=\$(samtools view -c -F 0x900 "\${BAM_FILE}")
  echo "QC | STATS | Total reads: \${TOTAL_READS}"

  # Mapped reads (primary, not unmapped)
  echo "QC | STATS | Counting mapped reads..."
  MAPPED_READS=\$(samtools view -c -F 0x904 "\${BAM_FILE}")
  echo "QC | STATS | Mapped reads: \${MAPPED_READS}"

  # Duplicate reads
  echo "QC | STATS | Counting duplicate reads..."
  DUP_READS=\$(samtools view -c -f 0x400 -F 0x900 "\${BAM_FILE}")
  echo "QC | STATS | Duplicate reads: \${DUP_READS}"

  # MAPQ filtered reads (primary, mapped, MAPQ≥threshold)
  echo "QC | STATS | Counting MAPQ≥\${MAPQ_THRESHOLD} reads..."
  MAPQ_READS=\$(samtools view -c -F 0x904 -q \${MAPQ_THRESHOLD} "\${BAM_FILE}")
  echo "QC | STATS | MAPQ≥\${MAPQ_THRESHOLD} reads: \${MAPQ_READS}"

  # MAPQ filtered + deduplicated reads
  echo "QC | STATS | Counting MAPQ≥\${MAPQ_THRESHOLD} + non-duplicate reads..."
  MAPQ_NODUP_READS=\$(samtools view -c -F 0xD04 -q \${MAPQ_THRESHOLD} "\${BAM_FILE}")
  echo "QC | STATS | MAPQ≥\${MAPQ_THRESHOLD} (no dup): \${MAPQ_NODUP_READS}"

  STATS_END=\$(date +%s)
  STATS_TIME=\$((STATS_END - STATS_START))
  echo "QC | STATS | Alignment statistics calculated in \${STATS_TIME}s"

  ###########################################################################
  # 5) STRAND BIAS ANALYSIS
  ###########################################################################

  echo "QC | STRAND | Analyzing strand bias..."

  STRAND_START=\$(date +%s)

  # Count reads on each strand (after MAPQ filtering)
  # Flag 0x10 = reverse strand
  samtools view -F 0x904 -q \${MAPQ_THRESHOLD} ${dedup_flag} "\${BAM_FILE}" | \
    awk '{
      if (and(\$2, 16)) {
        strand = "-"
      } else {
        strand = "+"
      }
      count[strand]++
    }
    END {
      print "strand\\tcount"
      print "+\\t" (count["+"] ? count["+"] : 0)
      print "-\\t" (count["-"] ? count["-"] : 0)
    }' > qc_strand_bias.tsv

  # Parse strand counts
  PLUS_READS=\$(awk 'NR==2 {print \$2}' qc_strand_bias.tsv)
  MINUS_READS=\$(awk 'NR==3 {print \$2}' qc_strand_bias.tsv)
  TOTAL_STRAND=\$((PLUS_READS + MINUS_READS))

  echo "QC | STRAND | Plus strand: \${PLUS_READS}"
  echo "QC | STRAND | Minus strand: \${MINUS_READS}"
  echo "QC | STRAND | Total (filtered): \${TOTAL_STRAND}"

  # Calculate strand balance
  if [[ \${TOTAL_STRAND} -gt 0 ]]; then
    PLUS_PERCENT=\$(awk -v p=\${PLUS_READS} -v t=\${TOTAL_STRAND} \
                       'BEGIN{printf "%.2f", 100.0*p/t}')
    MINUS_PERCENT=\$(awk -v m=\${MINUS_READS} -v t=\${TOTAL_STRAND} \
                        'BEGIN{printf "%.2f", 100.0*m/t}')
    echo "QC | STRAND | Plus: \${PLUS_PERCENT}%"
    echo "QC | STRAND | Minus: \${MINUS_PERCENT}%"
    
    # Check for severe bias (>70% or <30%)
    PLUS_INT=\${PLUS_PERCENT%.*}
    if [[ \${PLUS_INT} -gt 70 || \${PLUS_INT} -lt 30 ]]; then
      echo "QC | WARNING | Significant strand bias detected (expected ~50/50)"
    fi
  fi

  STRAND_END=\$(date +%s)
  STRAND_TIME=\$((STRAND_END - STRAND_START))
  echo "QC | STRAND | Strand analysis completed in \${STRAND_TIME}s"

  ###########################################################################
  # 6) FRAGMENT LENGTH DISTRIBUTION (PE only)
  ###########################################################################

  echo "QC | FRAGMENT | Analyzing fragment lengths..."

  FRAG_START=\$(date +%s)

  if [[ \${IS_PAIRED} -eq 1 ]]; then
    echo "QC | FRAGMENT | Extracting insert size distribution..."
    
    # Use samtools stats to get insert size distribution
    samtools stats -F 0x904 -q \${MAPQ_THRESHOLD} ${dedup_flag} "\${BAM_FILE}" | \
      awk '/^IS/ {print \$2 "\\t" \$3}' > frag_tmp.tsv || true
    
    if [[ -s frag_tmp.tsv ]]; then
      echo -e "fragment_length\\tcount" > qc_fragment_length.tsv
      cat frag_tmp.tsv >> qc_fragment_length.tsv
      rm -f frag_tmp.tsv
      
      FRAG_COUNT=\$(tail -n +2 qc_fragment_length.tsv | wc -l | tr -d ' ')
      echo "QC | FRAGMENT | Insert sizes: \${FRAG_COUNT} bins"
      
      # Calculate median and mean if possible
      MEDIAN_FRAG=\$(tail -n +2 qc_fragment_length.tsv | \
                     awk '{for(i=0;i<\$2;i++)print \$1}' | \
                     sort -n | \
                     awk '{a[NR]=\$1} END{print (NR%2==1)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}' || echo "NA")
      
      if [[ "\${MEDIAN_FRAG}" != "NA" ]]; then
        echo "QC | FRAGMENT | Median insert size: \${MEDIAN_FRAG} bp"
      fi
    else
      echo "QC | FRAGMENT | No fragment length data available"
      echo -e "fragment_length\\tcount" > qc_fragment_length.tsv
    fi
  else
    echo "QC | FRAGMENT | Single-end data, skipping fragment analysis"
    echo -e "fragment_length\\tcount" > qc_fragment_length.tsv
  fi

  FRAG_END=\$(date +%s)
  FRAG_TIME=\$((FRAG_END - FRAG_START))
  echo "QC | FRAGMENT | Fragment analysis completed in \${FRAG_TIME}s"

  ###########################################################################
  # 7) COVERAGE STATISTICS
  ###########################################################################

  echo "QC | COVERAGE | Calculating coverage statistics..."

  COV_START=\$(date +%s)

  # Use samtools coverage for fast per-chromosome coverage
  samtools coverage "\${BAM_FILE}" | \
    awk 'NR>1 {
      sum_cov += \$7
      sum_depth += \$7 * (\$3 - \$2)
      sum_len += (\$3 - \$2)
    }
    END {
      mean_depth = (sum_len > 0 ? sum_depth / sum_len : 0)
      print "metric\\tvalue"
      print "mean_coverage_depth\\t" mean_depth
    }' > qc_coverage.tsv

  MEAN_DEPTH=\$(awk 'NR==2 {print \$2}' qc_coverage.tsv)
  echo "QC | COVERAGE | Mean depth: \${MEAN_DEPTH}×"

  # Assess coverage quality
  DEPTH_INT=\${MEAN_DEPTH%.*}
  if [[ \${DEPTH_INT} -lt 10 ]]; then
    echo "QC | WARNING | Low coverage detected (<10×)"
  elif [[ \${DEPTH_INT} -ge 30 ]]; then
    echo "QC | COVERAGE | Good coverage (≥30×)"
  fi

  COV_END=\$(date +%s)
  COV_TIME=\$((COV_END - COV_START))
  echo "QC | COVERAGE | Coverage statistics calculated in \${COV_TIME}s"

  ###########################################################################
  # 8) PARSE UMI DEDUPLICATION STATISTICS
  ###########################################################################

  echo "QC | UMI | Checking UMI deduplication statistics..."

  UMI_ENABLED="false"
  UMI_INPUT_READS=0
  UMI_OUTPUT_READS=0
  UMI_DUPLICATES_REMOVED=0
  UMI_DEDUP_PERCENT=0

  if [[ -s "\${DEDUP_STATS}" ]]; then
    echo "QC | UMI | Parsing deduplication statistics..."
    
    # Priority 1: Our appended Summary section (reads_before=, reads_after= from samtools)
    # Use head -1 to avoid picking up stray matches; grep returns first match only
    if grep -q "=== Summary ===" "\${DEDUP_STATS}" 2>/dev/null; then
      UMI_ENABLED="true"
      UMI_INPUT_READS=\$(grep "^reads_before=" "\${DEDUP_STATS}" 2>/dev/null | head -1 | cut -d= -f2 || echo 0)
      UMI_OUTPUT_READS=\$(grep "^reads_after=" "\${DEDUP_STATS}" 2>/dev/null | head -1 | cut -d= -f2 || echo 0)
      UMI_DUPLICATES_REMOVED=\$(grep "^reads_removed=" "\${DEDUP_STATS}" 2>/dev/null | head -1 | cut -d= -f2 || echo 0)
      UMI_DEDUP_PERCENT=\$(grep "^percent_removed=" "\${DEDUP_STATS}" 2>/dev/null | head -1 | cut -d= -f2 || echo 0)
      
      # Sanity: output must not exceed input (indicates parse error or N/A)
      if [[ "\${UMI_INPUT_READS}" =~ ^[0-9]+\$ ]] && [[ "\${UMI_OUTPUT_READS}" =~ ^[0-9]+\$ ]] && [[ "\${UMI_OUTPUT_READS}" -gt "\${UMI_INPUT_READS}" ]] 2>/dev/null; then
        echo "QC | UMI | WARNING: reads_after > reads_before — parsing may have failed, trying umi_tools native format"
        UMI_INPUT_READS=0
        UMI_OUTPUT_READS=0
        UMI_DUPLICATES_REMOVED=0
        UMI_DEDUP_PERCENT=0
      fi
      # Also reset if we got N/A or non-numeric
      if [[ ! "\${UMI_INPUT_READS}" =~ ^[0-9]+\$ ]]; then
        UMI_INPUT_READS=0
      fi
      if [[ ! "\${UMI_OUTPUT_READS}" =~ ^[0-9]+\$ ]]; then
        UMI_OUTPUT_READS=0
      fi
      
      # If percentage not in file, calculate it
      if [[ "\${UMI_DEDUP_PERCENT}" == "0" || "\${UMI_DEDUP_PERCENT}" == "0.00" ]] && [[ \${UMI_INPUT_READS} -gt 0 ]] 2>/dev/null; then
        UMI_DEDUP_PERCENT=\$(awk -v d=\${UMI_DUPLICATES_REMOVED} \
                                 -v t=\${UMI_INPUT_READS} \
                                 'BEGIN{printf "%.2f", 100.0*d/t}')
      fi
    fi
    
    # Priority 2: Fallback to umi_tools native log format (Input Reads: or Input reads:, Number of reads out:)
    if [[ \${UMI_INPUT_READS:-0} -eq 0 || \${UMI_OUTPUT_READS:-0} -eq 0 ]] && grep -qiE "Input [Rr]eads:|Number of reads out:" "\${DEDUP_STATS}" 2>/dev/null; then
      UMI_ENABLED="true"
      UMI_INPUT_READS=\$(grep -i "Input Reads:" "\${DEDUP_STATS}" 2>/dev/null | tail -1 | sed -n 's/.*[Ii]nput [Rr]eads: *\\([0-9][0-9,]*\\).*/\\1/p' | tr -d ',' || echo 0)
      UMI_OUTPUT_READS=\$(grep -i "Number of reads out:" "\${DEDUP_STATS}" 2>/dev/null | tail -1 | sed -n 's/.*[Nn]umber of reads out: *\\([0-9][0-9,]*\\).*/\\1/p' | tr -d ',' || echo 0)
      if [[ \${UMI_INPUT_READS} -gt 0 ]] 2>/dev/null; then
        UMI_DUPLICATES_REMOVED=\$((UMI_INPUT_READS - UMI_OUTPUT_READS))
        UMI_DEDUP_PERCENT=\$(awk -v d=\${UMI_DUPLICATES_REMOVED} -v t=\${UMI_INPUT_READS} 'BEGIN{printf "%.2f", 100.0*d/t}')
      fi
      echo "QC | UMI | Parsed umi_tools native format (fallback)"
    fi
    
    if [[ "\${UMI_ENABLED}" == "true" ]]; then
      echo "QC | UMI | UMI deduplication enabled"
      echo "QC | UMI | Input reads: \${UMI_INPUT_READS}"
      echo "QC | UMI | Output reads: \${UMI_OUTPUT_READS}"
      echo "QC | UMI | Duplicates removed: \${UMI_DUPLICATES_REMOVED} (\${UMI_DEDUP_PERCENT}%)"
    else
      echo "QC | UMI | UMI deduplication not used"
    fi
  else
    echo "QC | UMI | No deduplication statistics available"
  fi

  ###########################################################################
  # 9) CALCULATE SUMMARY METRICS
  ###########################################################################

  echo "QC | METRICS | Calculating summary metrics..."

  # Calculate percentages
  if [[ \${TOTAL_READS} -gt 0 ]]; then
    MAP_PERCENT=\$(awk -v m=\${MAPPED_READS} -v t=\${TOTAL_READS} \
                      'BEGIN{printf "%.2f", 100.0*m/t}')
    DUP_PERCENT=\$(awk -v d=\${DUP_READS} -v t=\${TOTAL_READS} \
                      'BEGIN{printf "%.2f", 100.0*d/t}')
    MAPQ_PERCENT=\$(awk -v m=\${MAPQ_READS} -v t=\${TOTAL_READS} \
                       'BEGIN{printf "%.2f", 100.0*m/t}')
  else
    MAP_PERCENT=0
    DUP_PERCENT=0
    MAPQ_PERCENT=0
  fi

  # Calculate strand fraction
  if [[ \${TOTAL_STRAND} -gt 0 ]]; then
    PLUS_FRAC=\$(awk -v p=\${PLUS_READS} -v t=\${TOTAL_STRAND} \
                    'BEGIN{printf "%.4f", 1.0*p/t}')
  else
    PLUS_FRAC=0.5
  fi

  echo "QC | METRICS | Mapping rate: \${MAP_PERCENT}%"
  echo "QC | METRICS | Duplicate rate: \${DUP_PERCENT}%"
  echo "QC | METRICS | MAPQ≥\${MAPQ_THRESHOLD} pass rate: \${MAPQ_PERCENT}%"
  echo "QC | METRICS | Strand balance: \${PLUS_FRAC} (+) / \$(awk -v f=\${PLUS_FRAC} 'BEGIN{printf "%.4f", 1-f}') (-)"

  ###########################################################################
  # 10) WRITE JSON SUMMARY
  ###########################################################################

  echo "QC | OUTPUT | Writing JSON summary..."

  cat > qc_pol.json <<JSONEOF
{
  "sample_id": "\${SAMPLE_ID}",
  "condition": "\${CONDITION}",
  "timepoint": "\${TIMEPOINT}",
  "replicate": "\${REPLICATE}",
  "sequencing_mode": "\$([ \${IS_PAIRED} -eq 1 ] && echo "paired-end" || echo "single-end")",
  "total_reads_raw": \${TOTAL_READS},
  "mapped_reads": \${MAPPED_READS},
  "map_rate_percent": \${MAP_PERCENT},
  "duplicate_reads": \${DUP_READS},
  "duplicate_perc_of_total": \${DUP_PERCENT},
  "mapq_ge_\${MAPQ_THRESHOLD}_reads": \${MAPQ_READS},
  "mapq_ge_reads_nodup": \${MAPQ_NODUP_READS},
  "mapq_pass_percent": \${MAPQ_PERCENT},
  "strand_plus_reads": \${PLUS_READS},
  "strand_minus_reads": \${MINUS_READS},
  "strand_plus_fraction": \${PLUS_FRAC},
  "mean_coverage_depth": \${MEAN_DEPTH},
  "mapq_threshold": \${MAPQ_THRESHOLD},
  "deduplication_enabled": \${DEDUP_ENABLED},
  "umi_deduplication_enabled": \${UMI_ENABLED},
  "umi_input_reads": \${UMI_INPUT_READS},
  "umi_output_reads": \${UMI_OUTPUT_READS},
  "umi_duplicates_removed": \${UMI_DUPLICATES_REMOVED},
  "umi_deduplication_percent": \${UMI_DEDUP_PERCENT}\$([ \${IS_PAIRED} -eq 1 ] && echo ", \"median_fragment_length\": \${MEDIAN_FRAG:-null}" || echo "")
}
JSONEOF

  JSON_SIZE=\$(stat -c%s qc_pol.json 2>/dev/null || stat -f%z qc_pol.json 2>/dev/null || echo "unknown")
  echo "QC | OUTPUT | JSON summary: \${JSON_SIZE} bytes"

  ###########################################################################
  # 11) CREATE README
  ###########################################################################

  echo "QC | README | Creating documentation..."

  cat > README_qc.txt <<'DOCEOF'
================================================================================
QUALITY CONTROL REPORT — \${SAMPLE_ID}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Comprehensive quality control analysis for PRO-seq aligned reads.
  
  Sample Information:
    • Sample ID:  \${SAMPLE_ID}
    • Condition:  \${CONDITION}
    • Timepoint:  \${TIMEPOINT}
    • Replicate:  \${REPLICATE}
    • Mode:       \$([ \${IS_PAIRED} -eq 1 ] && echo "Paired-end" || echo "Single-end")

FILTERING CRITERIA
────────────────────────────────────────────────────────────────────────────
  Flags Excluded:
    • 0x100 (secondary alignment)
    • 0x800 (supplementary alignment)
    • 0x4   (unmapped)
    \$([ "\${DEDUP_ENABLED}" == "true" ] && echo "  • 0x400 (duplicate)" || echo "")
  
  Quality Threshold:
    • MAPQ ≥ \${MAPQ_THRESHOLD}

ALIGNMENT STATISTICS
────────────────────────────────────────────────────────────────────────────
  Total Reads:              \${TOTAL_READS}
  Mapped Reads:             \${MAPPED_READS} (\${MAP_PERCENT}%)
  Duplicate Reads:          \${DUP_READS} (\${DUP_PERCENT}%)
  MAPQ≥\${MAPQ_THRESHOLD} Reads:         \${MAPQ_READS} (\${MAPQ_PERCENT}%)
  MAPQ≥\${MAPQ_THRESHOLD} (no dup):      \${MAPQ_NODUP_READS}
  
  Quality Assessment:
    \$(if awk -v m=\${MAP_PERCENT%.*} 'BEGIN{exit (m>=70)?0:1}'; then
      echo "✓ Good mapping rate (≥70%)"
    elif awk -v m=\${MAP_PERCENT%.*} 'BEGIN{exit (m>=50)?0:1}'; then
      echo "⚠ Moderate mapping rate (50-70%)"
    else
      echo "✗ Low mapping rate (<50%)"
    fi)
    
    \$(if awk -v d=\${DUP_PERCENT%.*} 'BEGIN{exit (d<=30)?0:1}'; then
      echo "✓ Acceptable duplicate rate (≤30%)"
    else
      echo "⚠ High duplicate rate (>30%)"
    fi)

STRAND BIAS ANALYSIS
────────────────────────────────────────────────────────────────────────────
  Plus Strand (+):          \${PLUS_READS} (\${PLUS_PERCENT}%)
  Minus Strand (−):         \${MINUS_READS} (\${MINUS_PERCENT}%)
  
  Expected Balance:         ~50/50 (biology may cause slight imbalance)
  
  Quality Assessment:
    \$(if awk -v p=\${PLUS_PERCENT%.*} 'BEGIN{exit (p>=40 && p<=60)?0:1}'; then
      echo "✓ Good strand balance (40-60%)"
    elif awk -v p=\${PLUS_PERCENT%.*} 'BEGIN{exit (p>=30 && p<=70)?0:1}'; then
      echo "⚠ Moderate strand bias (30-70%)"
    else
      echo "✗ Severe strand bias (<30% or >70%)"
    fi)
  
  Note: PRO-seq libraries should show roughly equal strand distribution.
        Severe bias may indicate technical issues.

COVERAGE STATISTICS
────────────────────────────────────────────────────────────────────────────
  Mean Depth:               \${MEAN_DEPTH}×
  
  Expected Depth:           >10× (adequate), >30× (good)
  
  Quality Assessment:
    \$(if awk -v d=\${MEAN_DEPTH%.*} 'BEGIN{exit (d>=30)?0:1}'; then
      echo "✓ Good coverage (≥30×)"
    elif awk -v d=\${MEAN_DEPTH%.*} 'BEGIN{exit (d>=10)?0:1}'; then
      echo "✓ Adequate coverage (10-30×)"
    else
      echo "✗ Low coverage (<10×)"
    fi)

\$([ \${IS_PAIRED} -eq 1 ] && cat <<FRAGEOF
FRAGMENT LENGTH DISTRIBUTION (Paired-End)
────────────────────────────────────────────────────────────────────────────
  Median Insert Size:       \${MEDIAN_FRAG:-NA} bp
  
  Expected Range:           150-250 bp (nucleosome-sized)
  
  Quality Assessment:
    \$(if [[ "\${MEDIAN_FRAG}" != "NA" ]]; then
      MEDIAN_INT=\${MEDIAN_FRAG%.*}
      if [[ \${MEDIAN_INT} -ge 150 && \${MEDIAN_INT} -le 250 ]]; then
        echo "✓ Fragment size in expected range"
      else
        echo "⚠ Fragment size outside typical range"
      fi
    else
      echo "− Not available"
    fi)
  
  Note: PRO-seq typically shows enrichment for ~200 bp fragments
        corresponding to nucleosome-protected DNA.
FRAGEOF
)

\$([ "\${UMI_ENABLED}" == "true" ] && cat <<UMIEOF
UMI DEDUPLICATION STATISTICS
────────────────────────────────────────────────────────────────────────────
  UMI Status:               Enabled
  
  Input Reads:              \${UMI_INPUT_READS}
  Output Reads:             \${UMI_OUTPUT_READS}
  Duplicates Removed:       \${UMI_DUPLICATES_REMOVED} (\${UMI_DEDUP_PERCENT}%)
  
  Quality Assessment:
    \$(if awk -v d=\${UMI_DEDUP_PERCENT%.*} 'BEGIN{exit (d<=50)?0:1}'; then
      echo "✓ Reasonable duplication rate (≤50%)"
    elif awk -v d=\${UMI_DEDUP_PERCENT%.*} 'BEGIN{exit (d<=70)?0:1}'; then
      echo "⚠ High duplication rate (50-70%)"
    else
      echo "✗ Very high duplication rate (>70%)"
    fi)
  
  Note: UMI-based deduplication is more accurate than coordinate-based
        deduplication and can distinguish true duplicates from PCR copies.
UMIEOF
)

FILES GENERATED
────────────────────────────────────────────────────────────────────────────
  qc_strand_bias.tsv:
    Two-column TSV with strand (+/−) and read counts
    Used for: Assessing strand balance
  
  qc_fragment_length.tsv:
    Two-column TSV with fragment length and count
    Used for: Insert size distribution (PE only)
    Empty for single-end data
  
  qc_coverage.tsv:
    Two-column TSV with metric and value
    Currently contains: mean_coverage_depth
  
  qc_pol.json:
    JSON summary with all QC metrics
    Used for: Downstream aggregation and reporting
    Fields: 20+ metrics including alignment, strand, coverage stats
  
  README_qc.txt:
    This documentation file
  
  qc.log:
    Processing log with timestamps

INTERPRETATION GUIDE
────────────────────────────────────────────────────────────────────────────

Good Quality Sample:
  ✓ Mapping rate >70%
  ✓ Duplicate rate <30%
  ✓ Strand balance 40-60%
  ✓ Mean depth >10×
  ✓ Fragment size 150-250 bp (PE)

Potential Issues:
  ⚠ Low mapping rate (<70%):
    → Check adapter contamination
    → Verify genome reference
    → Check read quality
  
  ⚠ High duplicate rate (>30%):
    → May indicate low complexity library
    → Consider UMI-based deduplication
    → Check PCR cycles
  
  ⚠ Strand bias (<40% or >60%):
    → Check library preparation protocol
    → Verify PRO-seq specific steps
    → May indicate degradation
  
  ⚠ Low coverage (<10×):
    → Increase sequencing depth
    → Check library complexity
    → Verify enrichment efficiency

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These QC metrics are used for:
  
  1. Sample QC:
     - Identify failed libraries
     - Detect technical issues
     - Compare across batches
  
  2. Filtering:
     - Exclude low-quality samples
     - Identify outliers
     - Batch effect assessment
  
  3. Normalization:
     - Inform normalization strategy
     - Account for depth differences
     - Consider duplicate rates

PROCESSING DETAILS
────────────────────────────────────────────────────────────────────────────
  Tools Used:
    • samtools: Read counting, filtering, coverage
    • awk: Text processing and calculations
  
  Processing Times:
    • Alignment stats: \${STATS_TIME}s
    • Strand analysis: \${STRAND_TIME}s
    • Fragment analysis: \${FRAG_TIME}s
    • Coverage stats: \${COV_TIME}s
    • Total: \$((STATS_TIME + STRAND_TIME + FRAG_TIME + COV_TIME))s
  
  Samtools Flags Used:
    -F 0x900: Exclude secondary (0x100) and supplementary (0x800)
    -F 0x904: Above + exclude unmapped (0x4)
    -F 0xD04: Above + exclude duplicate (0x400) and unmapped (0x4)
    -q ${mapq_thr}: Minimum MAPQ threshold

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: \$(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: \${SAMPLE_ID}
  Module: 13_qc_pol_tracktx

================================================================================
DOCEOF

  echo "QC | README | Documentation created"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  TOTAL_TIME=\$((STATS_TIME + STRAND_TIME + FRAG_TIME + COV_TIME))

  echo "────────────────────────────────────────────────────────────────────────"
  echo "QC | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "QC | SUMMARY | Mode: \$([ \${IS_PAIRED} -eq 1 ] && echo "PE" || echo "SE")"
  echo "QC | SUMMARY | Total reads: \${TOTAL_READS}"
  echo "QC | SUMMARY | Mapped: \${MAPPED_READS} (\${MAP_PERCENT}%)"
  echo "QC | SUMMARY | Duplicates: \${DUP_READS} (\${DUP_PERCENT}%)"
  echo "QC | SUMMARY | MAPQ≥\${MAPQ_THRESHOLD}: \${MAPQ_READS} (\${MAPQ_PERCENT}%)"
  echo "QC | SUMMARY | Strand balance: \${PLUS_FRAC}"
  echo "QC | SUMMARY | Mean depth: \${MEAN_DEPTH}×"
  if [[ "\${UMI_ENABLED}" == "true" ]]; then
    echo "QC | SUMMARY | UMI deduplication: \${UMI_DEDUP_PERCENT}%"
  fi
  echo "QC | SUMMARY | Processing time: \${TOTAL_TIME}s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "QC | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}