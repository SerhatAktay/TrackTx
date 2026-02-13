// ============================================================================
// generate_tracks.nf — Strand-Specific Coverage Track Generation
// ============================================================================
//
// Purpose:
//   Generates nucleotide-resolution coverage tracks for PRO-seq analysis
//
// Critical Implementation Detail:
//   Uses bedtools genomecov -ibam for CORRECT 3'/5' end extraction
//   (bamToBed pipeline is INCORRECT as -3/-5 flags require BAM CIGAR data)
//
// Features:
//   • Direct BAM processing (no BED intermediate)
//   • Proper 3'/5' end extraction using alignment information
//   • Strand-specific tracks (positive and negative)
//   • Negative strand mirrored with -scale -1
//   • Both main BAM and allMap BAM processing
//   • Optional UMI deduplication
//   • bedGraph and BigWig output formats
//
// Coverage Types:
//   3' end: Always generated (PRO-seq standard)
//   5' end: Paired-end only (R2 alignment start positions)
//
// Inputs:
//   tuple(sample_id, filtered_bam, spikein_bam, condition, timepoint, replicate)
//   path(genome_fa)      : Genome FASTA (for chromosome sizes)
//   val(is_paired)       : "true" or "false"
//   path(allmap_bam)     : All mapped alignments BAM
//
// Outputs:
//   ${params.output_dir}/03_genome_tracks/${sample_id}/
//     ├── 3p/
//     │   ├── ${sample_id}.3p.pos.bedgraph       — Main BAM 3' positive
//     │   ├── ${sample_id}.3p.neg.bedgraph       — Main BAM 3' negative (mirrored)
//     │   ├── ${sample_id}.3p.pos.bw             — BigWig format
//     │   ├── ${sample_id}.3p.neg.bw             — BigWig format
//     │   ├── ${sample_id}.allMap.3p.pos.bedgraph
//     │   ├── ${sample_id}.allMap.3p.neg.bedgraph
//     │   ├── ${sample_id}.allMap.3p.pos.bw
//     │   └── ${sample_id}.allMap.3p.neg.bw
//     ├── 5p/ (PE only; SE creates empty placeholders)
//     │   └── [same structure as 3p/]
//     ├── ${sample_id}.dedup_stats.txt
//     ├── ${sample_id}.README_tracks.txt
//     └── tracks.log
//
// Parameters:
//   params.umi.enabled              : Enable UMI deduplication
//   params.umi.length               : UMI length in bases
//
// ============================================================================

nextflow.enable.dsl = 2

process generate_tracks {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/03_genome_tracks/${sample_id}",
             mode: 'copy',
             overwrite: true,
             saveAs: { filename ->
               def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
               // Exclude BAM files from publishing - they already exist in 02_alignments/
               // This prevents duplicating large BAM files in the tracks folder
               if (name.endsWith('.bam') || name.endsWith('.bam.bai')) {
                 return null  // Skip BAM files
               }
               return name
             }

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), path(filtered_bam), path(spikein_bam),
          val(condition), val(timepoint), val(replicate)
    path genome_fa
    val  is_paired
    path allmap_bam

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    // Main 3' tracks
    tuple val(sample_id),
          path("3p/${sample_id}.3p.pos.bedgraph"),
          path("3p/${sample_id}.3p.neg.bedgraph"),
          path("3p/${sample_id}.3p.pos.bw"),
          path("3p/${sample_id}.3p.neg.bw"),
          val(condition), val(timepoint), val(replicate),
          emit: bw3p_pair

    // Main 5' tracks (PE only)
    tuple val(sample_id),
          path("5p/${sample_id}.5p.pos.bedgraph", optional: true),
          path("5p/${sample_id}.5p.neg.bedgraph", optional: true),
          path("5p/${sample_id}.5p.pos.bw", optional: true),
          path("5p/${sample_id}.5p.neg.bw", optional: true),
          val(condition), val(timepoint), val(replicate),
          emit: bw5p_pair

    // allMap 3' tracks
    tuple val(sample_id),
          path("3p/${sample_id}.allMap.3p.pos.bedgraph"),
          path("3p/${sample_id}.allMap.3p.neg.bedgraph"),
          path("3p/${sample_id}.allMap.3p.pos.bw"),
          path("3p/${sample_id}.allMap.3p.neg.bw"),
          val(condition), val(timepoint), val(replicate),
          emit: allmap3p_pair

    // allMap 5' tracks (PE only)
    tuple val(sample_id),
          path("5p/${sample_id}.allMap.5p.pos.bedgraph", optional: true),
          path("5p/${sample_id}.allMap.5p.neg.bedgraph", optional: true),
          path("5p/${sample_id}.allMap.5p.pos.bw", optional: true),
          path("5p/${sample_id}.allMap.5p.neg.bw", optional: true),
          val(condition), val(timepoint), val(replicate),
          emit: allmap5p_pair

    // Legacy tuple for backward compatibility
    tuple val(sample_id),
          path(filtered_bam), path(spikein_bam),
          path("3p/${sample_id}.3p.pos.bedgraph"),
          path("3p/${sample_id}.3p.neg.bedgraph"),
          path("5p/${sample_id}.5p.pos.bedgraph", optional: true),
          path("5p/${sample_id}.5p.neg.bedgraph", optional: true),
          val(condition), val(timepoint), val(replicate),
          emit: track_tuple

    // BAM used for tracks (deduped when UMI on; same as density/pausing source)
    tuple val(sample_id), path("bam_for_downstream.bam"),
          val(condition), val(timepoint), val(replicate), emit: bam_for_tracks

    // Logs and documentation
    tuple val(sample_id), path("${sample_id}.dedup_stats.txt"),
          val(condition), val(timepoint), val(replicate), emit: dedup_stats
    path "${sample_id}.README_tracks.txt"
    path "tracks.log", emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a tracks.log) 2>&1

  # Error trap
  trap 'echo "TRACKS | ERROR | Process failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2; exit 1' ERR

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "TRACKS | START | sample=!{sample_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  echo "TRACKS | CONFIG | Initializing parameters..."

  SAMPLE_ID="!{sample_id}"
  THREADS=!{task.cpus}
  MAIN_BAM="!{filtered_bam}"
  ALLMAP_BAM="!{allmap_bam}"
  IS_PE="!{is_paired}"
  
  UMI_ENABLED="!{params.umi?.enabled ? 'true' : 'false'}"
  UMI_LENGTH=!{params.umi?.length ?: 0}

  # Constants
  readonly FLAG_PRIMARY_MAPPED=260      # Exclude unmapped(4) + secondary(256)
  readonly BIGWIG_TIMEOUT=600           # Seconds for BigWig conversion

  echo "TRACKS | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "TRACKS | CONFIG | Library type: $([ "${IS_PE}" == "true" ] && echo "Paired-end" || echo "Single-end")"
  echo "TRACKS | CONFIG | Threads: ${THREADS}"
  echo "TRACKS | CONFIG | UMI deduplication: ${UMI_ENABLED}"
  if [[ "${UMI_ENABLED}" == "true" ]]; then
    echo "TRACKS | CONFIG | UMI length: ${UMI_LENGTH} bp"
  fi

  ###########################################################################
  # 2) VALIDATION
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | VALIDATE | Checking inputs and environment..."
  echo "────────────────────────────────────────────────────────────────────────"

  VALIDATION_OK=1

  # Validate input BAMs
  if [[ ! -s "${MAIN_BAM}" ]]; then
    echo "TRACKS | ERROR | Main BAM file missing or empty: ${MAIN_BAM}"
    VALIDATION_OK=0
  else
    MAIN_SIZE=$(stat -c%s "${MAIN_BAM}" 2>/dev/null || stat -f%z "${MAIN_BAM}" 2>/dev/null || echo "unknown")
    echo "TRACKS | VALIDATE | Main BAM: ${MAIN_SIZE} bytes"
  fi

  if [[ ! -s "${ALLMAP_BAM}" ]]; then
    echo "TRACKS | ERROR | AllMap BAM file missing or empty: ${ALLMAP_BAM}"
    VALIDATION_OK=0
  else
    ALLMAP_SIZE=$(stat -c%s "${ALLMAP_BAM}" 2>/dev/null || stat -f%z "${ALLMAP_BAM}" 2>/dev/null || echo "unknown")
    echo "TRACKS | VALIDATE | AllMap BAM: ${ALLMAP_SIZE} bytes"
  fi

  # Validate required tools
  for TOOL in samtools bedtools bedGraphToBigWig; do
    if command -v ${TOOL} >/dev/null 2>&1; then
      echo "TRACKS | VALIDATE | ${TOOL}: $(command -v ${TOOL})"
    else
      echo "TRACKS | ERROR | Required tool not found: ${TOOL}"
      VALIDATION_OK=0
    fi
  done

  # Check for umi_tools if UMI dedup is enabled
  if [[ "${UMI_ENABLED}" == "true" ]]; then
    if command -v umi_tools >/dev/null 2>&1; then
      echo "TRACKS | VALIDATE | umi_tools: $(command -v umi_tools)"
    else
      echo "TRACKS | WARNING | umi_tools not found (will skip deduplication)"
    fi
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "TRACKS | ERROR | Validation failed"
    exit 1
  fi

  # Validate allMap BAM index (create if missing)
  if [[ ! -s "${ALLMAP_BAM}.bai" ]]; then
    echo "TRACKS | VALIDATE | Creating allMap BAM index..."
    samtools index -@ ${THREADS} "${ALLMAP_BAM}"
  fi

  echo "TRACKS | VALIDATE | All checks passed"

  # Report read counts
  MAIN_READS=$(samtools view -c -F ${FLAG_PRIMARY_MAPPED} "${MAIN_BAM}")
  ALLMAP_READS=$(samtools view -c -F ${FLAG_PRIMARY_MAPPED} "${ALLMAP_BAM}")
  echo "TRACKS | VALIDATE | Main BAM reads: ${MAIN_READS}"
  echo "TRACKS | VALIDATE | AllMap BAM reads: ${ALLMAP_READS}"

  # Create output directories
  mkdir -p 3p
  if [[ "${IS_PE}" == "true" ]]; then
    mkdir -p 5p
  fi
  echo "TRACKS | VALIDATE | Output directories created"

  ###########################################################################
  # 3) HELPER FUNCTIONS
  ###########################################################################

  # Convert bedGraph to BigWig with validation
  make_bigwig() {
    local bedgraph="$1"
    local bigwig="$2"
    
    if [[ ! -s "${bedgraph}" ]]; then
      echo "TRACKS | WARNING | Empty bedGraph, skipping BigWig: ${bedgraph}"
      return 0
    fi
    
    echo "TRACKS | BIGWIG | Converting: $(basename ${bedgraph}) → $(basename ${bigwig})"
    
    # Sort bedGraph (required for bedGraphToBigWig)
    echo "TRACKS | BIGWIG | Sorting bedGraph..."
    LC_ALL=C sort -k1,1 -k2,2n -o "${bedgraph}" "${bedgraph}"
    
    local line_count=$(wc -l < "${bedgraph}" | tr -d ' ')
    echo "TRACKS | BIGWIG | Sorted bedGraph: ${line_count} lines"
    
    # Convert to BigWig with timeout
    if ! timeout ${BIGWIG_TIMEOUT} bedGraphToBigWig "${bedgraph}" genome.sizes "${bigwig}"; then
      echo "TRACKS | ERROR | BigWig conversion failed or timed out: ${bedgraph}"
      return 1
    fi
    
    local bw_size=$(stat -c%s "${bigwig}" 2>/dev/null || stat -f%z "${bigwig}" 2>/dev/null || echo "unknown")
    echo "TRACKS | BIGWIG | Created: $(basename ${bigwig}) (${bw_size} bytes)"
    return 0
  }

  # Generate strand-specific coverage for specified end (3 or 5)
  # CRITICAL: Uses bedtools genomecov -ibam for correct 3'/5' end extraction
  # The -3/-5 flags require BAM CIGAR information (lost in BED conversion)
  generate_coverage() {
    local bam="$1"
    local end_type="$2"     # "3" or "5"
    local prefix="$3"       # Output file prefix
    
    echo "TRACKS | COVERAGE | Processing ${end_type}' end coverage for: $(basename ${bam})"
    
    # Validate BAM
    if ! samtools quickcheck "${bam}" 2>/dev/null; then
      echo "TRACKS | ERROR | Invalid BAM file: ${bam}"
      return 1
    fi
    
    # Positive strand coverage using direct BAM processing
    echo "TRACKS | COVERAGE | Computing positive strand..."
    if ! bedtools genomecov \
      -ibam "${bam}" \
      -${end_type} \
      -strand + \
      -bg \
      > "${prefix}.pos.bedgraph"; then
      echo "TRACKS | ERROR | Failed to generate positive strand coverage"
      return 1
    fi
    
    local pos_lines=$(wc -l < "${prefix}.pos.bedgraph" | tr -d ' ')
    local pos_size=$(stat -c%s "${prefix}.pos.bedgraph" 2>/dev/null || stat -f%z "${prefix}.pos.bedgraph" 2>/dev/null || echo "unknown")
    echo "TRACKS | COVERAGE | Positive strand: ${pos_lines} regions (${pos_size} bytes)"
    
    # Negative strand coverage (mirrored with -scale -1)
    echo "TRACKS | COVERAGE | Computing negative strand (mirrored with -scale -1)..."
    if ! bedtools genomecov \
      -ibam "${bam}" \
      -${end_type} \
      -strand - \
      -bg \
      -scale -1 \
      > "${prefix}.neg.bedgraph"; then
      echo "TRACKS | ERROR | Failed to generate negative strand coverage"
      return 1
    fi
    
    local neg_lines=$(wc -l < "${prefix}.neg.bedgraph" | tr -d ' ')
    local neg_size=$(stat -c%s "${prefix}.neg.bedgraph" 2>/dev/null || stat -f%z "${prefix}.neg.bedgraph" 2>/dev/null || echo "unknown")
    echo "TRACKS | COVERAGE | Negative strand: ${neg_lines} regions (${neg_size} bytes)"
    
    # Convert to BigWig
    echo "TRACKS | COVERAGE | Converting to BigWig format..."
    if ! make_bigwig "${prefix}.pos.bedgraph" "${prefix}.pos.bw"; then
      echo "TRACKS | ERROR | Failed to create positive strand BigWig"
      return 1
    fi
    if ! make_bigwig "${prefix}.neg.bedgraph" "${prefix}.neg.bw"; then
      echo "TRACKS | ERROR | Failed to create negative strand BigWig"
      return 1
    fi
    
    echo "TRACKS | COVERAGE | Complete: ${prefix}"
    return 0
  }

  # Perform UMI deduplication if enabled
  perform_umi_dedup() {
    local input_bam="$1"
    local output_bam="$2"
    local is_pe="$3"
    
    echo "TRACKS | DEDUP | Starting UMI deduplication..."
    echo "TRACKS | DEDUP | Mode: $([ "${is_pe}" == "true" ] && echo "Paired-end" || echo "Single-end")"
    
    if [[ "${is_pe}" == "true" ]]; then
      if ! umi_tools dedup \
        --paired \
        -I "${input_bam}" \
        -S "${output_bam}" \
        --log="${SAMPLE_ID}.dedup_stats.txt"; then
        echo "TRACKS | ERROR | umi_tools deduplication failed"
        return 1
      fi
    else
      if ! umi_tools dedup \
        -I "${input_bam}" \
        -S "${output_bam}" \
        --log="${SAMPLE_ID}.dedup_stats.txt"; then
        echo "TRACKS | ERROR | umi_tools deduplication failed"
        return 1
      fi
    fi
    
    # Index deduplicated BAM
    samtools index -@ ${THREADS} "${output_bam}"
    
    # Report statistics
    local before_reads=$(samtools view -c -F ${FLAG_PRIMARY_MAPPED} "${input_bam}")
    local after_reads=$(samtools view -c -F ${FLAG_PRIMARY_MAPPED} "${output_bam}")
    local removed=$((before_reads - after_reads))
    local pct_removed=$(awk -v b="${before_reads}" -v r="${removed}" 'BEGIN{printf "%.2f", (b>0)?(r*100.0/b):0}')
    
    echo "TRACKS | DEDUP | Reads before: ${before_reads}"
    echo "TRACKS | DEDUP | Reads after: ${after_reads}"
    echo "TRACKS | DEDUP | Removed: ${removed} (${pct_removed}%)"
    
    # Append summary to stats file
    {
      echo ""
      echo "=== Summary ==="
      echo "reads_before=${before_reads}"
      echo "reads_after=${after_reads}"
      echo "reads_removed=${removed}"
      echo "percent_removed=${pct_removed}"
    } >> "${SAMPLE_ID}.dedup_stats.txt"
    
    return 0
  }

  ###########################################################################
  # 4) OPTIONAL UMI DEDUPLICATION
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | DEDUP | Checking UMI deduplication settings..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Start with main BAM
  INPUT_BAM="${MAIN_BAM}"

  if [[ "${UMI_ENABLED}" == "true" && ${UMI_LENGTH} -gt 0 ]]; then
    if command -v umi_tools >/dev/null 2>&1; then
      # Create working copy and index
      cp "${MAIN_BAM}" aligned.bam
      samtools index -@ ${THREADS} aligned.bam
      
      # Perform deduplication
      if perform_umi_dedup "aligned.bam" "deduplicated.bam" "${IS_PE}"; then
        INPUT_BAM="deduplicated.bam"
        echo "TRACKS | DEDUP | Using deduplicated BAM for track generation"
      else
        echo "TRACKS | WARNING | Deduplication failed, using original BAM"
        {
          echo "umi_deduplication=failed"
          echo "status=umi_tools execution error"
          echo "reads_before=N/A"
          echo "reads_after=N/A"
          echo "reads_removed=N/A"
        } > "${SAMPLE_ID}.dedup_stats.txt"
      fi
    else
      echo "TRACKS | WARNING | umi_tools not found, skipping deduplication"
      {
        echo "umi_deduplication=skipped"
        echo "status=umi_tools not available"
        echo "reads_before=N/A"
        echo "reads_after=N/A"
        echo "reads_removed=N/A"
      } > "${SAMPLE_ID}.dedup_stats.txt"
    fi
  else
    echo "TRACKS | DEDUP | UMI deduplication disabled"
    {
      echo "umi_deduplication=disabled"
      echo "status=not requested in parameters"
      echo "reads_before=N/A"
      echo "reads_after=N/A"
      echo "reads_removed=N/A"
    } > "${SAMPLE_ID}.dedup_stats.txt"
  fi

  # Copy INPUT_BAM to named output for downstream (pol2 uses same BAM as tracks)
  echo "TRACKS | OUTPUT | Copying BAM used for tracks (deduped when UMI on)..."
  cp "${INPUT_BAM}" bam_for_downstream.bam
  samtools index -@ ${THREADS} bam_for_downstream.bam

  ###########################################################################
  # 5) PREPARE GENOME SIZES
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | GENOME | Extracting chromosome sizes..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Get chromosome sizes from BAM header (more reliable than FASTA)
  samtools view -H "${INPUT_BAM}" | \
    grep '^@SQ' | \
    cut -f2,3 | \
    sed 's/SN://g' | \
    sed 's/LN://g' \
    > genome.sizes

  if [[ ! -s genome.sizes ]]; then
    echo "TRACKS | ERROR | Failed to extract chromosome sizes from BAM header"
    exit 1
  fi

  CHR_COUNT=$(wc -l < genome.sizes | tr -d ' ')
  TOTAL_SIZE=$(awk '{sum+=$2} END{print sum}' genome.sizes)
  echo "TRACKS | GENOME | Chromosomes: ${CHR_COUNT}"
  echo "TRACKS | GENOME | Total genome size: ${TOTAL_SIZE} bp"

  # Display first few chromosomes
  echo "TRACKS | GENOME | First chromosomes:"
  head -5 genome.sizes | sed 's/^/TRACKS | GENOME |   /'

  ###########################################################################
  # 6) GENERATE 3' END COVERAGE (Always)
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | 3P | Generating 3' end coverage tracks..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Main BAM
  echo "TRACKS | 3P | Processing main BAM..."
  if ! generate_coverage "${INPUT_BAM}" "3" "3p/${SAMPLE_ID}.3p"; then
    echo "TRACKS | ERROR | Failed to generate 3' coverage from main BAM"
    exit 1
  fi

  # AllMap BAM
  echo "TRACKS | 3P | Processing allMap BAM..."
  if ! generate_coverage "${ALLMAP_BAM}" "3" "3p/${SAMPLE_ID}.allMap.3p"; then
    echo "TRACKS | ERROR | Failed to generate 3' coverage from allMap BAM"
    exit 1
  fi

  echo "TRACKS | 3P | 3' end coverage complete"

  ###########################################################################
  # 7) GENERATE 5' END COVERAGE (Paired-end only)
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  if [[ "${IS_PE}" == "true" ]]; then
    echo "TRACKS | 5P | Generating 5' end coverage tracks..."
    echo "────────────────────────────────────────────────────────────────────────"
    
    # Main BAM
    echo "TRACKS | 5P | Processing main BAM..."
    if ! generate_coverage "${INPUT_BAM}" "5" "5p/${SAMPLE_ID}.5p"; then
      echo "TRACKS | ERROR | Failed to generate 5' coverage from main BAM"
      exit 1
    fi
    
    # AllMap BAM
    echo "TRACKS | 5P | Processing allMap BAM..."
    if ! generate_coverage "${ALLMAP_BAM}" "5" "5p/${SAMPLE_ID}.allMap.5p"; then
      echo "TRACKS | ERROR | Failed to generate 5' coverage from allMap BAM"
      exit 1
    fi
    
    echo "TRACKS | 5P | 5' end coverage complete"
  else
    echo "TRACKS | 5P | Skipping 5' end coverage (single-end mode)"
    echo "────────────────────────────────────────────────────────────────────────"
    echo "TRACKS | 5P | Creating empty placeholder files for downstream compatibility"
    
    # Create 5p directory if not already created
    mkdir -p 5p
    
    # Create empty placeholder files (required for downstream modules)
    # Note: optional=true in output declarations doesn't work as intended
    touch "5p/${SAMPLE_ID}.5p.pos.bedgraph"
    touch "5p/${SAMPLE_ID}.5p.neg.bedgraph"
    touch "5p/${SAMPLE_ID}.5p.pos.bw"
    touch "5p/${SAMPLE_ID}.5p.neg.bw"
    touch "5p/${SAMPLE_ID}.allMap.5p.pos.bedgraph"
    touch "5p/${SAMPLE_ID}.allMap.5p.neg.bedgraph"
    touch "5p/${SAMPLE_ID}.allMap.5p.pos.bw"
    touch "5p/${SAMPLE_ID}.allMap.5p.neg.bw"
    
    echo "TRACKS | 5P | Placeholder files created"
  fi

  ###########################################################################
  # 8) CREATE DOCUMENTATION
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | README | Creating documentation..."
  echo "────────────────────────────────────────────────────────────────────────"

  cat > ${SAMPLE_ID}.README_tracks.txt <<'DOCEOF'
================================================================================
COVERAGE TRACKS — !{sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Strand-specific coverage tracks at nucleotide resolution for PRO-seq analysis.
  
  Track Types:
    • 3' end coverage: Always generated (PRO-seq standard)
    • 5' end coverage: Generated for paired-end data only
  
  BAM Sources:
    • Main BAM: Primary alignments !{params.umi?.enabled ? 'with UMI deduplication' : '(duplicates retained)'}
    • AllMap BAM: All mapped reads (primary + secondary alignments)

CRITICAL IMPLEMENTATION DETAIL
────────────────────────────────────────────────────────────────────────────
  This module uses bedtools genomecov -ibam for CORRECT 3'/5' end extraction.
  
  Why not bamToBed pipeline?
    The -3 and -5 flags require BAM CIGAR string information to determine
    true alignment end positions. BED format lacks this data, so:
    
    ❌ INCORRECT: bamToBed | bedtools genomecov -i stdin -3
       (produces wrong positions because CIGAR data is lost)
    
    ✅ CORRECT:   bedtools genomecov -ibam input.bam -3
       (properly extracts 3' positions from alignment data)

FILES
────────────────────────────────────────────────────────────────────────────

3' End Coverage (Always Generated):
  3p/!{sample_id}.3p.pos.bedgraph    — Positive strand (main BAM)
  3p/!{sample_id}.3p.neg.bedgraph    — Negative strand (main BAM, mirrored)
  3p/!{sample_id}.3p.pos.bw          — BigWig format (positive)
  3p/!{sample_id}.3p.neg.bw          — BigWig format (negative)
  
  3p/!{sample_id}.allMap.3p.pos.bedgraph — AllMap BAM positive strand
  3p/!{sample_id}.allMap.3p.neg.bedgraph — AllMap BAM negative strand (mirrored)
  3p/!{sample_id}.allMap.3p.pos.bw       — BigWig format
  3p/!{sample_id}.allMap.3p.neg.bw       — BigWig format

5' End Coverage (Paired-End Only):
  5p/!{sample_id}.5p.pos.bedgraph        — Main BAM positive strand
  5p/!{sample_id}.5p.neg.bedgraph        — Main BAM negative strand (mirrored)
  5p/!{sample_id}.5p.pos.bw              — BigWig format
  5p/!{sample_id}.5p.neg.bw              — BigWig format
  
  5p/!{sample_id}.allMap.5p.pos.bedgraph — AllMap BAM positive strand
  5p/!{sample_id}.allMap.5p.neg.bedgraph — AllMap BAM negative strand (mirrored)
  5p/!{sample_id}.allMap.5p.pos.bw       — BigWig format
  5p/!{sample_id}.allMap.5p.neg.bw       — BigWig format

  Note: For single-end data, 5p/ directory contains empty placeholder files
        for downstream pipeline compatibility.

Statistics:
  !{sample_id}.dedup_stats.txt       — UMI deduplication statistics
  tracks.log                         — Complete processing log

PROCESSING DETAILS
────────────────────────────────────────────────────────────────────────────

Pipeline Steps:
  1. Optional UMI deduplication (if enabled)
  2. Extract chromosome sizes from BAM header
  3. Generate coverage with bedtools genomecov -ibam
  4. Use -3 or -5 flags for end-specific coverage
  5. Negative strand multiplied by -1 using -scale flag
  6. Convert bedGraph to BigWig format

Command Example:
  # Positive strand 3' coverage
  bedtools genomecov -ibam sample.bam -3 -strand + -bg > pos.bedgraph
  
  # Negative strand 3' coverage (mirrored)
  bedtools genomecov -ibam sample.bam -3 -strand - -bg -scale -1 > neg.bedgraph

Key Settings:
  • Nucleotide resolution (single-base precision)
  • Strand-specific (separate positive and negative tracks)
  • Direct BAM processing (preserves alignment information)
  • Negative strand mirrored for UCSC Genome Browser compatibility
  • Chromosome sizes from BAM header (more reliable than FASTA)

UMI Deduplication:
  Status: !{params.umi?.enabled ? 'Enabled' : 'Disabled'}
  !{params.umi?.enabled ? 'Length: ' + params.umi.length + ' bp' : ''}
  !{params.umi?.enabled ? 'Duplicates removed before track generation' : 'Duplicates retained in coverage'}

USAGE
────────────────────────────────────────────────────────────────────────────

For Genome Browsers:
  • Load BigWig (.bw) files for visualization
  • Positive tracks show forward strand signal
  • Negative tracks show reverse strand signal (values < 0)
  • Pre-mirrored for direct UCSC Genome Browser viewing

For Computational Analysis:
  • Use bedGraph files for downstream processing
  • Values are RAW counts (not normalized at this stage)
  • Normalization occurs in subsequent pipeline steps
  • bedGraph format: chr<TAB>start<TAB>end<TAB>coverage

Main vs AllMap BAM Tracks:
  Main BAM tracks:
    • Primary alignments only
    • Cleaner signal
    • Use for most analyses
    • Recommended for peak calling
  
  AllMap BAM tracks:
    • Includes secondary alignments (multimappers)
    • Higher background signal
    • Use for multimapper-aware analyses
    • Used by divergent transcription detection

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These raw tracks will be:
  1. Normalized to CPM and siCPM (next module)
  2. Used for divergent transcription detection (allMap 3' tracks)
  3. Used for functional region calling (main 3' tracks)
  4. Used for Pol-II metrics calculation (normalized versions)

FILE FORMAT DETAILS
────────────────────────────────────────────────────────────────────────────

bedGraph Format:
  chromosome<TAB>start<TAB>end<TAB>coverage
  
  • Zero-based, half-open intervals [start, end)
  • Coverage values are raw read counts
  • Negative strand has negative values for visualization
  • Contiguous regions with same coverage are merged

BigWig Format:
  • Binary indexed format (faster than bedGraph)
  • Recommended for genome browser visualization
  • Created from bedGraph using UCSC bedGraphToBigWig
  • Allows efficient random access to genomic regions

QUALITY CHECKS
────────────────────────────────────────────────────────────────────────────

Expected Output:
  • Paired-end: 16 bedGraph + 16 BigWig files (32 total)
    - 8 files for 3' end coverage (4 main + 4 allMap)
    - 8 files for 5' end coverage (4 main + 4 allMap)
  • Single-end: 8 bedGraph + 8 BigWig files (16 total)
    - 8 files for 3' end coverage only
    - 5' directory contains empty placeholder files

Troubleshooting:
  • Empty bedGraph: Check if BAM has mapped reads
  • BigWig conversion failure: bedGraph sorting issue (automatically handled)
  • Missing 5' tracks: Normal for single-end data (placeholders created)
  • Large file sizes: Expected for high-coverage samples
  • Empty 5' placeholders: Normal for SE mode (downstream compatibility)

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • All coverage values are raw counts (not normalized)
  • Negative strand uses -scale -1 (native bedtools feature)
  • Chromosome sizes from BAM header ensure coordinate consistency
  • BigWig creation uses 600-second timeout for large genomes
  • bedGraph sorting is mandatory (automatically performed)
  • Empty placeholders ensure downstream module compatibility

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  UMI deduplication:    !{params.umi?.enabled ? 'Enabled' : 'Disabled'}
  UMI length:           !{params.umi?.length ?: 'N/A'} bp
  Library type:         !{is_paired == "true" ? "Paired-end" : "Single-end"}
  CPU threads:          !{task.cpus}

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Module:   06_generate_tracks
  Date:     $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample:   !{sample_id}

================================================================================
DOCEOF

  echo "TRACKS | README | Documentation created"

  ###########################################################################
  # 9) VALIDATION AND SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | VALIDATE | Verifying output files..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Count output files
  BG_COUNT=$(find 3p 5p -name "*.bedgraph" -type f 2>/dev/null | wc -l | tr -d ' ')
  BW_COUNT=$(find 3p 5p -name "*.bw" -type f 2>/dev/null | wc -l | tr -d ' ')

  # Expected file counts
  if [[ "${IS_PE}" == "true" ]]; then
    EXPECTED_BG=16  # 8 for 3p + 8 for 5p
    EXPECTED_BW=16
  else
    EXPECTED_BG=16  # 8 for 3p + 8 placeholders for 5p
    EXPECTED_BW=16
  fi

  echo "TRACKS | VALIDATE | bedGraph files: ${BG_COUNT}/${EXPECTED_BG}"
  echo "TRACKS | VALIDATE | BigWig files: ${BW_COUNT}/${EXPECTED_BW}"

  # Check critical 3' files (must be non-empty)
  CRITICAL_OK=1
  for file in \
    "3p/${SAMPLE_ID}.3p.pos.bedgraph" \
    "3p/${SAMPLE_ID}.3p.neg.bedgraph" \
    "3p/${SAMPLE_ID}.3p.pos.bw" \
    "3p/${SAMPLE_ID}.3p.neg.bw" \
    "3p/${SAMPLE_ID}.allMap.3p.pos.bedgraph" \
    "3p/${SAMPLE_ID}.allMap.3p.neg.bedgraph" \
    "3p/${SAMPLE_ID}.allMap.3p.pos.bw" \
    "3p/${SAMPLE_ID}.allMap.3p.neg.bw"; do
    
    if [[ ! -s "${file}" ]]; then
      echo "TRACKS | ERROR | Missing or empty critical file: ${file}"
      CRITICAL_OK=0
    else
      FILE_SIZE=$(stat -c%s "${file}" 2>/dev/null || stat -f%z "${file}" 2>/dev/null || echo "unknown")
      echo "TRACKS | VALIDATE | ${file}: ${FILE_SIZE} bytes"
    fi
  done

  # Check that 5' files exist (may be empty for SE)
  for file in \
    "5p/${SAMPLE_ID}.5p.pos.bedgraph" \
    "5p/${SAMPLE_ID}.5p.neg.bedgraph" \
    "5p/${SAMPLE_ID}.5p.pos.bw" \
    "5p/${SAMPLE_ID}.5p.neg.bw" \
    "5p/${SAMPLE_ID}.allMap.5p.pos.bedgraph" \
    "5p/${SAMPLE_ID}.allMap.5p.neg.bedgraph" \
    "5p/${SAMPLE_ID}.allMap.5p.pos.bw" \
    "5p/${SAMPLE_ID}.allMap.5p.neg.bw"; do
    
    if [[ ! -f "${file}" ]]; then
      echo "TRACKS | ERROR | Missing file: ${file}"
      CRITICAL_OK=0
    fi
  done

  if [[ ${CRITICAL_OK} -eq 0 ]]; then
    echo "TRACKS | ERROR | Validation failed"
    exit 1
  fi

  echo "TRACKS | VALIDATE | All required files present"

  ###########################################################################
  # FINAL SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | SUMMARY | Processing Complete"
  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "TRACKS | SUMMARY | Library type: $([ "${IS_PE}" == "true" ] && echo "PE (3' + 5')" || echo "SE (3' only)")"
  echo "TRACKS | SUMMARY | UMI dedup: ${UMI_ENABLED}"
  echo "TRACKS | SUMMARY | bedGraph files: ${BG_COUNT}"
  echo "TRACKS | SUMMARY | BigWig files: ${BW_COUNT}"
  echo "TRACKS | SUMMARY | Total output size: $(du -sh . 2>/dev/null | cut -f1 || echo "unknown")"
  
  if [[ "${IS_PE}" != "true" ]]; then
    echo "TRACKS | SUMMARY | Note: 5' tracks are placeholders (SE mode)"
  fi
  
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "TRACKS | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}
