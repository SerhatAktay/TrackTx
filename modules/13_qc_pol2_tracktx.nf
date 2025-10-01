// ============================================================================
// qc_pol2_tracktx.nf — library QC for Pol II–type datasets (DSL2)
// ----------------------------------------------------------------------------
// Overview
//   • Filters: primary, mapped, pass-QC, (optional) non-dup, MAPQ ≥ threshold
//   • Strand bias (+ / −) from filtered reads
//   • Fragment-length histogram (PE only) via `samtools stats IS`
//   • Mean coverage using `samtools coverage` (fast) + idxstats for genome length
//   • JSON summary with totals, mapped, dups, MAPQ-pass %, strand fractions, mean depth
//
// Inputs
//   tuple( sample_id, aligned_bam, condition, timepoint, replicate )
//
// Outputs (publishDir):
//   ${params.output_dir}/10_qc/${sample_id}/
//     ├── qc_strand_bias.tsv
//     ├── qc_fragment_length.tsv
//     ├── qc_coverage.tsv
//     ├── qc_pol2.json
//     └── README.txt
// ============================================================================

nextflow.enable.dsl = 2

process qc_pol2_tracktx {

  // Meta / resources
  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/10_qc/${sample_id}", mode: 'copy', overwrite: true

  // Tooling
  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // Inputs
  input:
    tuple val(sample_id),
          path(aligned_bam),
          val(condition), val(timepoint), val(replicate)

  // Declared outputs
  output:
    tuple val(sample_id), path('qc_strand_bias.tsv'), path('qc_fragment_length.tsv'), path('qc_coverage.tsv'), val(condition), val(timepoint), val(replicate)
    path 'qc_pol2.json', emit: json
    tuple val(sample_id), path('qc_pol2.json'), val(condition), val(timepoint), val(replicate), emit: json_meta
    path 'README.txt'

  script:
  def mapq_thr = params.qc?.mapq ?: 10
  def dedup_flag = (params.qc?.dedup == null || params.qc.dedup) ? '-F 0x400' : ''
  """
  set -euo pipefail
  
  # ── 1) Strand bias (MAPQ filtered, primary, mapped, no dup if requested) ──
  samtools view -F 0x904 -q ${mapq_thr} ${dedup_flag} "${aligned_bam}" \
    | awk '{if(and(\$2,16)) strand="-"; else strand="+"; count[strand]++} 
           END{print "strand\\tcount"; 
               print "+\\t" (count["+"]?count["+"]:0); 
               print "-\\t" (count["-"]?count["-"]:0)}' \
    > qc_strand_bias.tsv

  # ── 2) Fragment length histogram (PE only) ────────────────────────────────
  if samtools view -H "${aligned_bam}" | grep -q '@PG.*bowtie2'; then
    # PE data
    samtools stats -F 0x904 -q ${mapq_thr} ${dedup_flag} "${aligned_bam}" \
      | awk '/^IS/ {print \$2 "\\t" \$3}' > frag_tmp.tsv || true
    if [[ -s frag_tmp.tsv ]]; then
      echo -e "fragment_length\\tcount" > qc_fragment_length.tsv
      cat frag_tmp.tsv >> qc_fragment_length.tsv
    else
      echo -e "fragment_length\\tcount" > qc_fragment_length.tsv
    fi
  else
    # SE data
    echo -e "fragment_length\\tcount" > qc_fragment_length.tsv
  fi

  # ── 3) Coverage stats ──────────────────────────────────────────────────────
  samtools coverage "${aligned_bam}" \
    | awk 'NR>1 {sum_cov+=\$7; sum_depth+=\$7*(\$3-\$2); sum_len+=(\$3-\$2)} 
           END{mean_depth = (sum_len>0 ? sum_depth/sum_len : 0); 
               print "metric\\tvalue"; 
               print "mean_coverage_depth\\t" mean_depth}' \
    > qc_coverage.tsv

  # ── 4) Build QC JSON with full metrics ─────────────────────────────────────
  # Total reads
  total_reads=\$(samtools view -c -F 0x900 "${aligned_bam}")
  
  # Mapped reads (primary, not unmapped)
  mapped_reads=\$(samtools view -c -F 0x904 "${aligned_bam}")
  
  # Duplicate reads
  dup_reads=\$(samtools view -c -f 0x400 -F 0x900 "${aligned_bam}")
  
  # MAPQ filtered reads (primary, mapped, MAPQ≥threshold)
  mapq_reads=\$(samtools view -c -F 0x904 -q ${mapq_thr} "${aligned_bam}")
  
  # MAPQ filtered + deduplicated reads
  mapq_nodup_reads=\$(samtools view -c -F 0xD04 -q ${mapq_thr} "${aligned_bam}")
  
  # Strand counts
  plus_reads=\$(awk 'NR==2 {print \$2}' qc_strand_bias.tsv)
  minus_reads=\$(awk 'NR==3 {print \$2}' qc_strand_bias.tsv)
  
  # Mean depth
  mean_depth=\$(awk 'NR==2 {print \$2}' qc_coverage.tsv)
  
  # Calculate percentages
  map_percent=\$(awk -v m=\$mapped_reads -v t=\$total_reads 'BEGIN{print (t>0 ? 100.0*m/t : 0)}')
  dup_percent=\$(awk -v d=\$dup_reads -v t=\$total_reads 'BEGIN{print (t>0 ? 100.0*d/t : 0)}')
  mapq_percent=\$(awk -v m=\$mapq_reads -v t=\$total_reads 'BEGIN{print (t>0 ? 100.0*m/t : 0)}')
  plus_frac=\$(awk -v p=\$plus_reads -v t=\$((\$plus_reads + \$minus_reads)) 'BEGIN{print (t>0 ? 1.0*p/t : 0.5)}')
  
  # Write JSON
  cat > qc_pol2.json <<JSON
{
  "sample_id": "${sample_id}",
  "total_reads_raw": \$total_reads,
  "mapped_reads": \$mapped_reads,
  "map_rate_percent": \$map_percent,
  "duplicate_reads": \$dup_reads,
  "duplicate_perc_of_total": \$dup_percent,
  "mapq_ge_${mapq_thr}_reads": \$mapq_reads,
  "mapq_ge_reads_nodup": \$mapq_nodup_reads,
  "dedup_reads_mapq_ge": \$mapq_nodup_reads,
  "mapq_pass_percent": \$mapq_percent,
  "strand_plus_reads": \$plus_reads,
  "strand_minus_reads": \$minus_reads,
  "strand_plus_fraction": \$plus_frac,
  "mean_coverage_depth": \$mean_depth,
  "mapq_threshold": ${mapq_thr},
  "deduplication_enabled": ${params.qc?.dedup ?: true}
}
JSON

  # ── 5) README ──────────────────────────────────────────────────────────────
  cat > README.txt <<TXT
TrackTx QC for sample: ${sample_id}

Filters Applied:
  • Primary alignments only (excludes secondary/supplementary)
  • Mapped reads only
  • MAPQ ≥ ${mapq_thr}
  • ${params.qc?.dedup ?: true ? 'Duplicates removed' : 'Duplicates included'}

Outputs:
  • qc_strand_bias.tsv       - Read counts by strand
  • qc_fragment_length.tsv   - Insert size distribution (PE only)
  • qc_coverage.tsv          - Mean coverage depth
  • qc_pol2.json             - Summary metrics in JSON format

Metrics:
  • Total reads: \$total_reads
  • Mapped: \$mapped_reads (\$map_percent%)
  • Duplicates: \$dup_reads (\$dup_percent%)
  • MAPQ≥${mapq_thr} (deduplicated): \$mapq_nodup_reads
  • Mean depth: \$mean_depth
TXT
  """
}
