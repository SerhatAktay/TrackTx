// ============================================================================
// run_alignment.nf — Bowtie2 Alignment with Spike-in Support
// ============================================================================
//
// Purpose:
//   Performs Bowtie2 alignment for PRO-seq data with proper strand handling
//   and optional spike-in normalization support
//
// Key Features:
//   • PRO-seq convention: R1 reverse-complemented before alignment
//   • Paired-end support with proper mate orientation (--ff)
//   • Primary genome alignment → filtered and unfiltered BAMs
//   • Optional spike-in alignment of unmapped reads (SE mode)
//   • Comprehensive QC metrics and summaries
//
// Inputs:
//   tuple(sample_id, read1, read2, condition, timepoint, replicate)
//   tuple(genome_id, source, genome_fa) + path genome_bt2
//   tuple(spike_id, spike_source, spike_fa) + path spike_bt2
//
// Outputs:
//   ${params.output_dir}/02_alignments/${sample_id}/
//     ├── ${sample_id}_allMap.bam        — All alignments (primary + secondary)
//     ├── ${sample_id}.bam/.bai          — Primary alignments only
//     ├── ${sample_id}_spikein.bam/.bai  — Spike-in alignments
//     ├── *.flagstat / *.idxstats
//     ├── bowtie2_*.log
//     ├── aligner_summary.tsv
//     ├── insert_size.tsv (PE only)
//     ├── README_alignment.txt
//     └── align_reads.log
//
// ============================================================================

nextflow.enable.dsl = 2

process run_alignment {

  // ── Process Configuration ────────────────────────────────────────────────
  tag        { sample_id }
  label      'conda'
  cache      'deep'  // Ignore resource allocation changes for better caching

  publishDir "${params.output_dir}/02_alignments/${sample_id}",
             mode: 'copy',
             overwrite: true

  // ── Inputs ───────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), path(read1), path(read2),
          val(condition), val(timepoint), val(replicate)

    tuple val(genome_id), val(source), path(genome_fa)
    path  genome_bt2

    tuple val(spike_id), val(spike_source), path(spike_fa)
    path  spike_bt2
    
    val(is_paired_end)  // Pass as input to avoid params hash pollution

  // ── Outputs ──────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("${sample_id}.bam"),
          path("${sample_id}_allMap.bam"),
          path("${sample_id}_spikein.bam"),
          val(condition), val(timepoint), val(replicate)

    path "bowtie2_*.log",               emit: align_logs
    path "*.flagstat",                  emit: flagstats
    path "*.idxstats",                  emit: idxstats
    path "aligner_summary.tsv",         emit: align_summary
    path "insert_size.tsv", optional: true, emit: insert_size

    path "README_alignment.txt"
    path "align_reads.log",             emit: log

  // ── Main Script ──────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to unified log
  exec > >(tee -a align_reads.log) 2>&1

  # Error trap
  trap 'echo "ALIGN | ERROR | Process failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2; exit 1' ERR

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "ALIGN | START | sample=!{sample_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  echo "ALIGN | CONFIG | Initializing parameters..."

  SAMPLE_ID="!{sample_id}"
  GENOME_ID="!{genome_id}"
  SPIKE_ID="!{spike_id}"
  
  R1="!{read1}"
  R2="!{read2}"

  THREADS=!{task.cpus}
  # Threading strategy: use all threads for Bowtie2, minimal for samtools sorting
  BT2_THREADS="${THREADS}"
  SAM_THREADS=1

  # Determine if paired-end from input val and R2 presence
  PAIRED_PARAM='!{is_paired_end ? "true" : "false"}'
  IS_PE="false"
  if [[ "${PAIRED_PARAM}" == "true" && -s "${R2}" ]]; then
    IS_PE="true"
  fi

  echo "ALIGN | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "ALIGN | CONFIG | Genome: ${GENOME_ID}"
  echo "ALIGN | CONFIG | Spike-in: ${SPIKE_ID}"
  echo "ALIGN | CONFIG | Library type: $([ "${IS_PE}" == "true" ] && echo "Paired-end" || echo "Single-end")"
  echo "ALIGN | CONFIG | Threads: ${THREADS} (Bowtie2: ${BT2_THREADS}, samtools: ${SAM_THREADS})"
  echo "ALIGN | CONFIG | R1: ${R1}"
  if [[ "${IS_PE}" == "true" ]]; then
    echo "ALIGN | CONFIG | R2: ${R2}"
  fi

  # Validate paired-end configuration
  if [[ "${PAIRED_PARAM}" == "true" && "${IS_PE}" != "true" ]]; then
    echo "ALIGN | ERROR | params.paired_end=true but R2 missing or empty" >&2
    exit 1
  fi

  ###########################################################################
  # 2) VALIDATION
  ###########################################################################

  echo "ALIGN | VALIDATE | Checking system resources and inputs..."

  # Check disk space
  AVAIL_GB=$(df -BG . 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G' || echo "0")
  echo "ALIGN | VALIDATE | Available disk space: ${AVAIL_GB} GB"
  if [[ "$AVAIL_GB" -lt 10 ]]; then
    echo "ALIGN | WARNING | Low disk space: ${AVAIL_GB}GB available"
    echo "ALIGN | WARNING | Alignment may fail with large datasets"
  fi

  # Check input files
  VALIDATION_OK=1

  for file in "${R1}"; do
    if [[ ! -s "${file}" ]]; then
      echo "ALIGN | ERROR | Read file missing or empty: ${file}"
      VALIDATION_OK=0
    else
      FILE_SIZE=$(stat -c%s "${file}" 2>/dev/null || stat -f%z "${file}" 2>/dev/null || echo "unknown")
      echo "ALIGN | VALIDATE | R1 size: ${FILE_SIZE} bytes"
    fi
  done

  if [[ "${IS_PE}" == "true" && ! -s "${R2}" ]]; then
    echo "ALIGN | ERROR | R2 file missing or empty: ${R2}"
    VALIDATION_OK=0
  elif [[ "${IS_PE}" == "true" ]]; then
    FILE_SIZE=$(stat -c%s "${R2}" 2>/dev/null || stat -f%z "${R2}" 2>/dev/null || echo "unknown")
    echo "ALIGN | VALIDATE | R2 size: ${FILE_SIZE} bytes"
  fi

  # Validate required tools
  for TOOL in bowtie2 samtools gzip; do
    if command -v ${TOOL} >/dev/null 2>&1; then
      echo "ALIGN | VALIDATE | ${TOOL}: $(command -v ${TOOL})"
    else
      echo "ALIGN | ERROR | Required tool not found: ${TOOL}"
      VALIDATION_OK=0
    fi
  done

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "ALIGN | ERROR | Validation failed"
    exit 1
  fi

  echo "ALIGN | VALIDATE | All checks passed"

  ###########################################################################
  # 3) HELPER FUNCTIONS
  ###########################################################################

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  else
    PYTHON_CMD="python3"
  fi

  # Decompress input files (handles .gz and uncompressed)
  decompress() {
    [[ "$1" == *.gz ]] && gzip -cd -- "$1" || cat -- "$1"
  }

  # Reverse-complement FASTQ stream (PRO-seq convention for R1)
  rc_stream() {
    if command -v seqkit >/dev/null 2>&1; then
      seqkit seq --quiet -t dna -r -p -j "${THREADS}"
    else
      ${PYTHON_CMD} - "$@" <<'PYEND'
import sys
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
it = iter(sys.stdin)
while True:
    try:
        h = next(it).rstrip('\n')
        s = next(it).rstrip('\n')
        p = next(it).rstrip('\n')
        q = next(it).rstrip('\n')
    except StopIteration:
        break
    sys.stdout.write(f"{h}\n{s.translate(comp)[::-1]}\n{p}\n{q[::-1]}\n")
PYEND
    fi
  }

  ###########################################################################
  # 4) STAGE BOWTIE2 INDICES
  ###########################################################################

  echo "ALIGN | INDEX | Staging Bowtie2 indices..."

  mkdir -p bt2_index

  # Stage genome index
  GENOME_INDEX_LIST='!{(genome_bt2 instanceof List) ? genome_bt2.collect{it.toString()}.join(' ') : genome_bt2.toString()}'
  if [[ -n "${GENOME_INDEX_LIST}" ]]; then
    # shellcheck disable=SC2086
    cp -f ${GENOME_INDEX_LIST} bt2_index/
    echo "ALIGN | INDEX | Genome index files staged"
  fi
  GENOME_IDX="bt2_index/${GENOME_ID}"

  # Stage spike-in index (if provided)
  SPIKE_INDEX_LIST='!{(spike_bt2 instanceof List) ? spike_bt2.collect{it.toString()}.join(' ') : spike_bt2.toString()}'
  SPIKE_IDX=""
  if [[ -n "${SPIKE_ID}" && "${SPIKE_ID}" != "none" && -n "${SPIKE_INDEX_LIST}" ]]; then
    # shellcheck disable=SC2086
    if ! cp -f ${SPIKE_INDEX_LIST} bt2_index/ 2>/dev/null; then
      echo "ALIGN | ERROR | Failed to stage spike-in index files" >&2
      echo "ALIGN | DEBUG | SPIKE_INDEX_LIST: ${SPIKE_INDEX_LIST}" >&2
      ls -la bt2_index/ 2>/dev/null || true
      exit 1
    fi
    SPIKE_IDX="bt2_index/${SPIKE_ID}"
    echo "ALIGN | INDEX | Spike-in index files staged"
  else
    echo "ALIGN | INDEX | No spike-in index provided"
  fi

  # Verify spike-in index completeness (when used)
  if [[ -n "${SPIKE_IDX}" ]]; then
    SPIKE_INDEX_COMPLETE=0
    for ext in bt2 bt2l; do
      ALL_PRESENT=1
      for shard in 1 2 3 4 rev.1 rev.2; do
        if [[ ! -s "${SPIKE_IDX}.${shard}.${ext}" ]]; then
          ALL_PRESENT=0
          break
        fi
      done
      if [[ ${ALL_PRESENT} -eq 1 ]]; then
        SPIKE_INDEX_COMPLETE=1
        echo "ALIGN | INDEX | Spike-in index verified (${ext})"
        break
      fi
    done
    if [[ ${SPIKE_INDEX_COMPLETE} -eq 0 ]]; then
      echo "ALIGN | ERROR | Incomplete spike-in index: ${SPIKE_IDX}" >&2
      echo "ALIGN | DEBUG | Files in bt2_index/:" >&2
      ls -lh bt2_index/ 2>/dev/null || true
      exit 1
    fi
  fi

  # Verify genome index completeness
  echo "ALIGN | INDEX | Verifying index integrity..."
  INDEX_COMPLETE=0
  for ext in bt2 bt2l; do
    ALL_PRESENT=1
    for shard in 1 2 3 4 rev.1 rev.2; do
      if [[ ! -s "${GENOME_IDX}.${shard}.${ext}" ]]; then
        ALL_PRESENT=0
        break
      fi
    done
    if [[ ${ALL_PRESENT} -eq 1 ]]; then
      INDEX_COMPLETE=1
      echo "ALIGN | INDEX | Found complete ${ext} index family"
      break
    fi
  done

  if [[ ${INDEX_COMPLETE} -eq 0 ]]; then
    echo "ALIGN | ERROR | No complete Bowtie2 index found for: ${GENOME_IDX}"
    echo "ALIGN | ERROR | Expected shards: .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2"
    echo "ALIGN | DEBUG | Files present in bt2_index/:"
    ls -lh bt2_index/ || true
    exit 1
  fi

  ###########################################################################
  # 5) PRIMARY GENOME ALIGNMENT
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | PRIMARY | Starting genome alignment (${GENOME_ID})..."
  echo "────────────────────────────────────────────────────────────────────────"

  if [[ "${IS_PE}" == "true" ]]; then
    echo "ALIGN | PRIMARY | Mode: Paired-end with PRO-seq orientation (--ff)"
    echo "ALIGN | PRIMARY | -1: original R2, -2: RC(R1)"
    
    bowtie2 -p "${BT2_THREADS}" \
            --end-to-end \
            --ff \
            --no-unal \
            -x "${GENOME_IDX}" \
            -1 <(decompress "${R2}") \
            -2 <(decompress "${R1}" | rc_stream) \
            --un-conc unaligned_R%.fastq \
            2> bowtie2_primary.log \
    | samtools sort -@ "${SAM_THREADS}" -o "${SAMPLE_ID}_allMap.bam"
    
    # Combine unaligned reads for spike-in
    cat unaligned_R1.fastq unaligned_R2.fastq > unaligned.fastq
    
  else
    echo "ALIGN | PRIMARY | Mode: Single-end with RC(R1)"
    
    bowtie2 -p "${BT2_THREADS}" \
            --end-to-end \
            --no-unal \
            -x "${GENOME_IDX}" \
            -U <(decompress "${R1}" | rc_stream) \
            --un unaligned.fastq \
            2> bowtie2_primary.log \
    | samtools sort -@ "${SAM_THREADS}" -o "${SAMPLE_ID}_allMap.bam"
  fi

  # Verify allMap BAM
  samtools quickcheck -v "${SAMPLE_ID}_allMap.bam"
  ALLMAP_SIZE=$(stat -c%s "${SAMPLE_ID}_allMap.bam" 2>/dev/null || stat -f%z "${SAMPLE_ID}_allMap.bam" 2>/dev/null || echo "unknown")
  echo "ALIGN | PRIMARY | allMap BAM created: ${ALLMAP_SIZE} bytes"

  ###########################################################################
  # 6) EXTRACT PRIMARY ALIGNMENTS
  ###########################################################################

  echo "ALIGN | FILTER | Extracting primary alignments only..."
  echo "ALIGN | FILTER | Filtering: -F 260 (unmapped=4, secondary=256)"

  # Extract primary mapped reads (exclude unmapped=4 and secondary=256)
  samtools view -@ "${THREADS}" -h -b -F 260 "${SAMPLE_ID}_allMap.bam" \
  | samtools sort -@ "${SAM_THREADS}" -o "${SAMPLE_ID}.bam"

  # Verify and index primary BAM
  samtools quickcheck -v "${SAMPLE_ID}.bam"
  samtools index -@ "${THREADS}" "${SAMPLE_ID}.bam"
  
  PRIMARY_SIZE=$(stat -c%s "${SAMPLE_ID}.bam" 2>/dev/null || stat -f%z "${SAMPLE_ID}.bam" 2>/dev/null || echo "unknown")
  echo "ALIGN | FILTER | Primary BAM created: ${PRIMARY_SIZE} bytes"
  echo "ALIGN | FILTER | BAM indexed successfully"

  # Generate QC metrics
  echo "ALIGN | QC | Generating flagstat and idxstats..."
  samtools flagstat -@ "${THREADS}" "${SAMPLE_ID}.bam" > "${SAMPLE_ID}.flagstat"
  samtools idxstats "${SAMPLE_ID}.bam" > "${SAMPLE_ID}.idxstats"

  ###########################################################################
  # 7) SPIKE-IN ALIGNMENT
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | SPIKEIN | Processing spike-in alignment..."
  echo "────────────────────────────────────────────────────────────────────────"

  if [[ -n "${SPIKE_IDX}" && -s unaligned.fastq ]]; then
    UNALIGNED_COUNT=$(wc -l < unaligned.fastq | awk '{print $1/4}')
    echo "ALIGN | SPIKEIN | Aligning ${UNALIGNED_COUNT} unaligned reads to ${SPIKE_ID}"
    echo "ALIGN | SPIKEIN | Mode: Single-end (all spike-in alignments as SE)"
    
    # Use -m 4G for samtools sort to avoid OOM with large unaligned sets (e.g. >20M reads)
    if ! bowtie2 -p "${BT2_THREADS}" \
            --end-to-end \
            --no-unal \
            -x "${SPIKE_IDX}" \
            -U unaligned.fastq \
            2> bowtie2_spikein.log \
    | samtools sort -@ "${SAM_THREADS}" -m 4G -o "${SAMPLE_ID}_spikein.bam"; then
      echo "ALIGN | SPIKEIN | ERROR | Spike-in alignment failed" >&2
      echo "ALIGN | SPIKEIN | DEBUG | bowtie2_spikein.log:" >&2
      cat bowtie2_spikein.log 2>/dev/null | tail -100 >&2
      exit 1
    fi
    
    SPIKEIN_SIZE=$(stat -c%s "${SAMPLE_ID}_spikein.bam" 2>/dev/null || stat -f%z "${SAMPLE_ID}_spikein.bam" 2>/dev/null || echo "unknown")
    echo "ALIGN | SPIKEIN | Spike-in BAM created: ${SPIKEIN_SIZE} bytes"
    
  else
    if [[ -z "${SPIKE_IDX}" ]]; then
      echo "ALIGN | SPIKEIN | No spike-in index provided - creating empty BAM"
    else
      echo "ALIGN | SPIKEIN | No unaligned reads - creating empty BAM"
    fi
    
    # Create empty but valid BAM with proper header
    samtools view -@ "${THREADS}" -H "${SAMPLE_ID}.bam" \
    | samtools view -@ "${THREADS}" -b - > "${SAMPLE_ID}_spikein.bam"
    : > bowtie2_spikein.log
  fi

  # Verify and index spike-in BAM
  samtools quickcheck -v "${SAMPLE_ID}_spikein.bam"
  samtools index -@ "${THREADS}" "${SAMPLE_ID}_spikein.bam"
  
  echo "ALIGN | SPIKEIN | Generating spike-in QC metrics..."
  samtools flagstat -@ "${THREADS}" "${SAMPLE_ID}_spikein.bam" > "${SAMPLE_ID}_spikein.flagstat"
  samtools idxstats "${SAMPLE_ID}_spikein.bam" > "${SAMPLE_ID}_spikein.idxstats"

  ###########################################################################
  # 8) GENERATE ALIGNMENT SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | SUMMARY | Generating alignment statistics..."
  echo "────────────────────────────────────────────────────────────────────────"

  {
    echo -e "metric\tvalue"
    
    # Primary genome alignment metrics
    tot=$(awk '/in total/ {print $1}' "${SAMPLE_ID}.flagstat")
    map=$(awk '/ mapped \\(/ {print $1}' "${SAMPLE_ID}.flagstat")
    mpr=$(awk -F'[()% ]+' '/ mapped \\(/ {print $(NF-1)}' "${SAMPLE_ID}.flagstat")
    sec=$(awk '/secondary/ {print $1}' "${SAMPLE_ID}.flagstat")
    dup=$(awk '/duplicates/ {print $1}' "${SAMPLE_ID}.flagstat")
    
    # Ensure numeric values
    [[ -n "$tot" ]] || tot=0
    [[ -n "$map" ]] || map=0
    [[ -n "$mpr" ]] || mpr=0
    [[ -n "$sec" ]] || sec=0
    [[ -n "$dup" ]] || dup=0
    
    echo -e "genome_total_reads\t${tot}"
    echo -e "genome_mapped_reads\t${map}"
    echo -e "genome_map_rate_pct\t${mpr}"
    echo -e "genome_secondary_reads\t${sec}"
    echo -e "genome_duplicate_reads\t${dup}"
    
    # Spike-in alignment metrics
    stot=$(awk '/in total/ {print $1}' "${SAMPLE_ID}_spikein.flagstat")
    smap=$(awk '/ mapped \\(/ {print $1}' "${SAMPLE_ID}_spikein.flagstat")
    smpr=$(awk -F'[()% ]+' '/ mapped \\(/ {print $(NF-1)}' "${SAMPLE_ID}_spikein.flagstat")
    
    [[ -n "$stot" ]] || stot=0
    [[ -n "$smap" ]] || smap=0
    [[ -n "$smpr" ]] || smpr=0
    
    echo -e "spike_total_reads\t${stot}"
    echo -e "spike_mapped_reads\t${smap}"
    echo -e "spike_map_rate_pct\t${smpr}"
    
    # Paired-end specific metrics (concordance)
    if [[ "${IS_PE}" == "true" && -s bowtie2_primary.log ]]; then
      conc0=$(awk '/aligned concordantly 0 times/ {print $1}' bowtie2_primary.log)
      conc1=$(awk '/aligned concordantly exactly 1 time/ {print $1}' bowtie2_primary.log)
      concM=$(awk '/aligned concordantly >1 times/ {print $1}' bowtie2_primary.log)
      dis=$(awk '/aligned discordantly 1 time/ {print $1}' bowtie2_primary.log)
      
      [[ -n "$conc0" ]] || conc0=0
      [[ -n "$conc1" ]] || conc1=0
      [[ -n "$concM" ]] || concM=0
      [[ -n "$dis" ]] || dis=0
      
      echo -e "pe_concordant_0\t${conc0}"
      echo -e "pe_concordant_1\t${conc1}"
      echo -e "pe_concordant_gt1\t${concM}"
      echo -e "pe_discordant_1\t${dis}"
    fi
    
  } > aligner_summary.tsv

  echo "ALIGN | SUMMARY | Alignment summary created"

  # Extract insert size distribution for paired-end
  if [[ "${IS_PE}" == "true" ]]; then
    echo "ALIGN | SUMMARY | Extracting insert size distribution..."
    samtools stats "${SAMPLE_ID}.bam" \
    | awk -F'\\t' '$1=="IS"{print $2"\\t"$3}' \
    > insert_size.tsv || true
  fi

  ###########################################################################
  # 9) GENERATE DOCUMENTATION
  ###########################################################################

  echo "ALIGN | README | Creating documentation..."

  cat > README_alignment.txt <<'DOCEOF'
================================================================================
ALIGNMENT ARTIFACTS — PRO-seq Pipeline
================================================================================

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample ID:     !{sample_id}
  Genome:        !{genome_id}
  Spike-in:      !{spike_id}
  Library type:  !{is_paired_end ? "Paired-end" : "Single-end"}

OUTPUT FILES
────────────────────────────────────────────────────────────────────────────

Primary Alignments:
  !{sample_id}_allMap.bam          — All mapped reads (primary + secondary)
  !{sample_id}.bam                 — Primary alignments only (-F 260)
  !{sample_id}.bam.bai             — BAM index

Spike-in Alignments:
  !{sample_id}_spikein.bam         — Spike-in alignments (SE mode)
  !{sample_id}_spikein.bam.bai     — BAM index

Quality Control:
  !{sample_id}.flagstat            — Primary BAM statistics
  !{sample_id}.idxstats            — Per-chromosome alignment counts
  !{sample_id}_spikein.flagstat    — Spike-in BAM statistics
  !{sample_id}_spikein.idxstats    — Spike-in per-chromosome counts

Alignment Logs:
  bowtie2_primary.log              — Bowtie2 primary alignment log
  bowtie2_spikein.log              — Bowtie2 spike-in alignment log

Summary Reports:
  aligner_summary.tsv              — Key alignment metrics (TSV format)
  insert_size.tsv                  — Insert size distribution (PE only)

Documentation:
  README_alignment.txt             — This file
  align_reads.log                  — Complete processing log (stdout + stderr)

PROCESSING DETAILS
────────────────────────────────────────────────────────────────────────────

PRO-seq Convention:
  • R1 is reverse-complemented before alignment
  • This captures the 3' end of nascent RNA at polymerase position
  • R2 (if PE) is aligned as-is

Paired-End Alignment:
  • Uses --ff orientation flag (both mates on forward strand)
  • Effective mate configuration: -1 original_R2, -2 RC(R1)
  • Unaligned read pairs combined for spike-in alignment

Single-End Alignment:
  • Aligns RC(R1) to genome
  • Unaligned reads saved for spike-in alignment

Spike-in Alignment:
  • Always performed in single-end mode
  • Uses reads that failed to map to primary genome
  • Empty BAM created if no spike-in index or no unmapped reads

Filtering Strategy:
  • sample_allMap.bam: Contains ALL alignments (primary + secondary)
  • sample.bam: Primary alignments only (excludes unmapped=4, secondary=256)
  • Duplicates are RETAINED at this stage (removed in deduplication step)

ALIGNMENT METRICS
────────────────────────────────────────────────────────────────────────────

aligner_summary.tsv contains:
  • genome_total_reads        : Total reads processed
  • genome_mapped_reads       : Successfully mapped reads
  • genome_map_rate_pct       : Mapping percentage
  • genome_secondary_reads    : Secondary/multimapper reads
  • genome_duplicate_reads    : PCR/optical duplicates
  • spike_total_reads         : Reads attempted for spike-in
  • spike_mapped_reads        : Spike-in mapped reads
  • spike_map_rate_pct        : Spike-in mapping percentage

Paired-end additional metrics:
  • pe_concordant_0          : Pairs aligned 0 times concordantly
  • pe_concordant_1          : Pairs aligned exactly 1 time concordantly
  • pe_concordant_gt1        : Pairs aligned >1 times concordantly
  • pe_discordant_1          : Pairs aligned discordantly

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────

These BAM files are used by subsequent pipeline modules:
  1. generate_tracks.nf     — Creates bedGraph/BigWig coverage tracks
  2. normalize_tracks.nf    — CPM and siCPM normalization
  3. detect_divergent.nf    — Divergent transcription detection
  4. call_regions.nf        — Functional region calling

BAM File Selection Guide:
  • Use sample.bam for most analyses (cleaner signal)
  • Use sample_allMap.bam for multimapper-aware analyses
  • Use sample_spikein.bam for spike-in normalization calculations

QUALITY CHECKS
────────────────────────────────────────────────────────────────────────────

Expected Metrics:
  • Mapping rate >70% for well-prepared libraries
  • Secondary reads <10% (genome-dependent)
  • Spike-in mapping rate variable (depends on spike-in abundance)

Troubleshooting:
  • Low mapping rate: Check library quality, adapter contamination
  • High secondary rate: May indicate repetitive regions or multimappers
  • Empty spike-in BAM: Normal if no spike-in index provided
  • Missing insert_size.tsv: Normal for single-end data

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────

Bowtie2 Parameters:
  • --end-to-end : Align entire read (no soft clipping)
  • --no-unal    : Suppress unaligned reads in output
  • --ff         : PE orientation (both mates forward strand)

SAM Flags:
  • 4   : Read unmapped
  • 256 : Secondary alignment (multimapper)
  • -F 260 : Exclude both unmapped and secondary reads

File Naming:
  • _allMap : Contains all mappings (primary + secondary)
  • No suffix: Primary alignments only
  • _spikein : Spike-in genome alignments

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  Genome index:       !{genome_id}
  Spike-in index:     !{spike_id}
  Paired-end mode:    !{is_paired_end}
  CPU threads:        !{task.cpus}

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Module:   05_run_alignment
  Date:     $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample:   !{sample_id}

================================================================================
DOCEOF

  echo "ALIGN | README | Documentation complete"

  ###########################################################################
  # 10) FINAL VALIDATION AND SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | VALIDATE | Verifying all output files..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Check critical output files
  FINAL_OK=1
  for file in \
    "${SAMPLE_ID}_allMap.bam" \
    "${SAMPLE_ID}.bam" \
    "${SAMPLE_ID}.bam.bai" \
    "${SAMPLE_ID}_spikein.bam" \
    "${SAMPLE_ID}_spikein.bam.bai" \
    "${SAMPLE_ID}.flagstat" \
    "${SAMPLE_ID}.idxstats" \
    "aligner_summary.tsv"; do
    
    if [[ ! -s "${file}" ]]; then
      echo "ALIGN | ERROR | Missing or empty output file: ${file}"
      FINAL_OK=0
    else
      SIZE=$(stat -c%s "${file}" 2>/dev/null || stat -f%z "${file}" 2>/dev/null || echo "unknown")
      echo "ALIGN | VALIDATE | ${file}: ${SIZE} bytes"
    fi
  done

  if [[ ${FINAL_OK} -eq 0 ]]; then
    echo "ALIGN | ERROR | Output validation failed"
    exit 1
  fi

  echo "ALIGN | VALIDATE | All outputs verified"

  # Print summary statistics
  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | SUMMARY | Alignment Statistics"
  echo "────────────────────────────────────────────────────────────────────────"
  
  if command -v column >/dev/null 2>&1; then
    column -t -s $'\\t' aligner_summary.tsv | sed 's/^/ALIGN | SUMMARY | /'
  else
    sed 's/^/ALIGN | SUMMARY | /' aligner_summary.tsv
  fi
  
  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | SUMMARY | Output Files"
  echo "────────────────────────────────────────────────────────────────────────"
  echo "ALIGN | SUMMARY | Total size: $(du -sh . 2>/dev/null | cut -f1 || echo "unknown")"
  echo "ALIGN | SUMMARY | BAM files: $(find . -name "*.bam" -type f | wc -l | tr -d ' ')"
  echo "ALIGN | SUMMARY | Index files: $(find . -name "*.bai" -type f | wc -l | tr -d ' ')"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "ALIGN | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}