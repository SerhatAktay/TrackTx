// ============================================================================
// generate_coverage_tracks.nf — Strand-Specific Coverage Track Generation
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
//   5' end: Always generated (PE and SE; bedtools genomecov -5)
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
//     ├── 5p/
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


process generate_coverage_tracks {

  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  // Persistent storage for cross-run caching
  storeDir   { "${params.assets_dir ?: "${projectDir}/assets"}/03_genome_tracks/${sample_id}" }

  publishDir { "${params.output_dir}/03_genome_tracks/${sample_id}" },
             mode: params.publish_mode,
             overwrite: true,
             saveAs: { filename ->
               def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
               def pathStr = filename.toString()
               // Exclude BAM files from publishing - they already exist in 02_alignments/
               if (name.endsWith('.bam') || name.endsWith('.bam.bai')) return null
               // Skip entire folder when output.raw_tracks: false (~150 MB saved)
               if (params.get('output')?.get('raw_tracks')?.toString() == 'false') return null
               // Skip 5' tracks when norm.emit_5p: false (~75 MB saved)
               if (params.get('norm')?.get('emit_5p')?.toString() == 'false' && (pathStr.contains('5p/') || name.contains('.5p.'))) return null
               // Skip allMap tracks when norm.emit_allmap: false (~25 MB in 03)
               if (params.get('norm')?.get('emit_allmap')?.toString() == 'false' && name.contains('allMap')) return null
               // Skip bedGraphs when output.bedgraph: false (BigWigs sufficient for genome browsers)
               if (params.get('output')?.get('bedgraph')?.toString() == 'false' && name.endsWith('.bedgraph')) return null
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

    // Main 5' tracks
    tuple val(sample_id),
          path("5p/${sample_id}.5p.pos.bedgraph"),
          path("5p/${sample_id}.5p.neg.bedgraph"),
          path("5p/${sample_id}.5p.pos.bw"),
          path("5p/${sample_id}.5p.neg.bw"),
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

    // allMap 5' tracks
    tuple val(sample_id),
          path("5p/${sample_id}.allMap.5p.pos.bedgraph"),
          path("5p/${sample_id}.allMap.5p.neg.bedgraph"),
          path("5p/${sample_id}.allMap.5p.pos.bw"),
          path("5p/${sample_id}.allMap.5p.neg.bw"),
          val(condition), val(timepoint), val(replicate),
          emit: allmap5p_pair

    // Legacy tuple for backward compatibility
    tuple val(sample_id),
          path(filtered_bam), path(spikein_bam),
          path("3p/${sample_id}.3p.pos.bedgraph"),
          path("3p/${sample_id}.3p.neg.bedgraph"),
          path("5p/${sample_id}.5p.pos.bedgraph"),
          path("5p/${sample_id}.5p.neg.bedgraph"),
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
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  # umi_tools imports matplotlib; avoid font-cache stall when UMI dedup is enabled
  export MPLCONFIGDIR="\${TMPDIR:-/tmp}/matplotlib"

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a tracks.log)
  exec 2> >(tee -a tracks.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "generate_coverage_tracks" "Unexpected process failure" "Check tracks.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "TRACKS | START | sample=${sample_id} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  echo "TRACKS | CONFIG | Initializing parameters..."

  SAMPLE_ID="${sample_id}"
  THREADS=${task.cpus}
  MAIN_BAM="${filtered_bam}"
  ALLMAP_BAM="${allmap_bam}"
  IS_PE="${is_paired}"
  
  UMI_ENABLED="${params.umi?.enabled ? 'true' : 'false'}"
  UMI_LENGTH=${params.umi?.length ?: 0}

  # Constants
  readonly FLAG_PRIMARY_MAPPED=260      # Exclude unmapped(4) + secondary(256)
  readonly BIGWIG_TIMEOUT=600           # Seconds for BigWig conversion

  # Cross-sample I/O lock: on this pipeline's typical single USB/HDD-backed
  # work volume, running several samples' large whole-BAM sequential reads
  # (index builds, full-file copies) at once doesn't parallelize disk
  # throughput -- it interleaves them into seek-heavy access. Measured live:
  # 3 concurrent 5GB allMap-BAM index builds ran at ~7MB/s EACH (12+min for
  # a file that should take ~1min at a 5400rpm 2.5" HDD's real ~80-130MB/s
  # single-stream ceiling). flock here serializes just those whole-file
  # passes across concurrently running samples/tasks; it does NOT touch the
  # CPU-bound work (chromosome-sharded dedup, genomecov|sort, bigwig) that
  # was already fixed to use its allotted cpus. Lock lives on the shared
  # work volume (projectDir, bind-mounted identically into every task's
  # container) so every concurrent task sees the same lock file. Override
  # slot-wait timeout with TRACKS_IO_LOCK_TIMEOUT; a timeout fails the task
  # through the normal ERR trap below rather than silently racing unlocked.
  IO_LOCK_FILE="${projectDir}/.tracktx_io.lock"
  IO_LOCK_TIMEOUT=\${TRACKS_IO_LOCK_TIMEOUT:-1800}
  with_io_lock() {
    flock -w "\${IO_LOCK_TIMEOUT}" "\${IO_LOCK_FILE}" "\$@"
  }

  echo "TRACKS | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "TRACKS | CONFIG | Library type: \$([ "\${IS_PE}" == "true" ] && echo "Paired-end" || echo "Single-end")"
  echo "TRACKS | CONFIG | Threads: \${THREADS}"
  echo "TRACKS | CONFIG | UMI deduplication: \${UMI_ENABLED}"
  if [[ "\${UMI_ENABLED}" == "true" ]]; then
    echo "TRACKS | CONFIG | UMI length: \${UMI_LENGTH} bp"
  fi

  ###########################################################################
  # 2) VALIDATION
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | VALIDATE | Checking inputs and environment..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Validate input BAMs
  if [[ ! -s "\${MAIN_BAM}" ]]; then
    tracktx_error "generate_coverage_tracks" "Main BAM missing or empty: \${MAIN_BAM}" "Check align_reads_to_genome produced sample.bam"
  fi
  MAIN_SIZE=\$(stat -c%s "\${MAIN_BAM}" 2>/dev/null || stat -f%z "\${MAIN_BAM}" 2>/dev/null || echo "unknown")
  echo "TRACKS | VALIDATE | Main BAM: \${MAIN_SIZE} bytes"

  if [[ ! -s "\${ALLMAP_BAM}" ]]; then
    tracktx_error "generate_coverage_tracks" "AllMap BAM missing or empty: \${ALLMAP_BAM}" "Check align_reads_to_genome produced sample_allMap.bam"
  fi
  ALLMAP_SIZE=\$(stat -c%s "\${ALLMAP_BAM}" 2>/dev/null || stat -f%z "\${ALLMAP_BAM}" 2>/dev/null || echo "unknown")
  echo "TRACKS | VALIDATE | AllMap BAM: \${ALLMAP_SIZE} bytes"

  # Validate required tools
  for TOOL in samtools bedtools bedGraphToBigWig; do
    if ! command -v \${TOOL} >/dev/null 2>&1; then
      tracktx_error "generate_coverage_tracks" "Required tool not found: \${TOOL}" "Install \${TOOL} or use -profile docker"
    fi
    echo "TRACKS | VALIDATE | \${TOOL}: \$(command -v \${TOOL})"
  done

  # Check for umi_tools if UMI dedup is enabled
  if [[ "\${UMI_ENABLED}" == "true" ]]; then
    if command -v umi_tools >/dev/null 2>&1; then
      echo "TRACKS | VALIDATE | umi_tools: \$(command -v umi_tools)"
    else
      echo "TRACKS | WARNING | umi_tools not found (will skip deduplication)"
    fi
  fi

  # Validate allMap BAM index (create if missing)
  if [[ ! -s "\${ALLMAP_BAM}.bai" ]]; then
    echo "TRACKS | VALIDATE | Creating allMap BAM index..."
    with_io_lock samtools index -@ \${THREADS} "\${ALLMAP_BAM}"
  fi

  echo "TRACKS | VALIDATE | All checks passed"

  # Report read counts
  MAIN_READS=\$(samtools view -c -F \${FLAG_PRIMARY_MAPPED} "\${MAIN_BAM}")
  ALLMAP_READS=\$(samtools view -c -F \${FLAG_PRIMARY_MAPPED} "\${ALLMAP_BAM}")
  echo "TRACKS | VALIDATE | Main BAM reads: \${MAIN_READS}"
  echo "TRACKS | VALIDATE | AllMap BAM reads: \${ALLMAP_READS}"

  # Create output directories
  mkdir -p 3p 5p
  echo "TRACKS | VALIDATE | Output directories created"

  ###########################################################################
  # 3) HELPER FUNCTIONS
  ###########################################################################

  # Convert bedGraph to BigWig with validation
  make_bigwig() {
    local bedgraph="\$1"
    local bigwig="\$2"
    
    if [[ ! -s "\${bedgraph}" ]]; then
      echo "TRACKS | WARNING | Empty bedGraph, skipping BigWig: \${bedgraph}"
      return 0
    fi
    
    echo "TRACKS | BIGWIG | Converting: \$(basename \${bedgraph}) → \$(basename \${bigwig})"

    # NOTE: bedGraphs reaching here are ALREADY coordinate-sorted — generate_coverage
    # pipes "bedtools genomecov | sort" so the on-disk bedGraph is written once,
    # already sorted. This removes a full write+read+rewrite of every (often
    # >1 GB) bedGraph on the work disk, which is the dominant cost when work/ is
    # on a slow/USB-backed volume. Do NOT re-sort here.
    local line_count=\$(wc -l < "\${bedgraph}" | tr -d ' ')
    echo "TRACKS | BIGWIG | Pre-sorted bedGraph: \${line_count} lines"

    # Convert to BigWig with timeout
    if ! timeout \${BIGWIG_TIMEOUT} bedGraphToBigWig "\${bedgraph}" genome.sizes "\${bigwig}"; then
      echo "TRACKS | ERROR | BigWig conversion failed or timed out: \${bedgraph}"
      return 1
    fi
    
    local bw_size=\$(stat -c%s "\${bigwig}" 2>/dev/null || stat -f%z "\${bigwig}" 2>/dev/null || echo "unknown")
    echo "TRACKS | BIGWIG | Created: \$(basename \${bigwig}) (\${bw_size} bytes)"
    return 0
  }

  # Generate strand-specific coverage for specified end (3 or 5)
  # CRITICAL: Uses bedtools genomecov -ibam for correct 3'/5' end extraction
  # The -3/-5 flags require BAM CIGAR information (lost in BED conversion)
  generate_coverage() {
    local bam="\$1"
    local end_type="\$2"     # "3" or "5"
    local prefix="\$3"       # Output file prefix
    
    echo "TRACKS | COVERAGE | Processing \${end_type}' end coverage for: \$(basename \${bam})"
    
    # Validate BAM
    if ! samtools quickcheck "\${bam}" 2>/dev/null; then
      echo "TRACKS | ERROR | Invalid BAM file: \${bam}"
      return 1
    fi
    
    # Sort budget: cap per-sort memory so it stays in RAM rather than spilling
    # thousands of tiny temp files (catastrophic on slow/USB work dirs), without
    # oversubscribing the task's RAM. SORT_MEM and the concurrency cap (MAX_PAR)
    # are computed together in the concurrency section below and exported, so this
    # default is just a standalone fallback. Temp goes to a fast dir (never the
    # USB-backed work dir). Override with SORT_MEM / SORT_TMPDIR.
    : "\${SORT_MEM:=\$(( ${task.memory.toGiga()} * 50 / 100 / 2 ))G}"
    : "\${SORT_TMP:=\${SORT_TMPDIR:-/tmp}}"
    : "\${SORT_PARALLEL:=\${THREADS}}"
    mkdir -p "\${SORT_TMP}" 2>/dev/null || SORT_TMP=/tmp

    # Positive strand coverage, sorted in a single streamed pass.
    # genomecov | sort writes the (already-sorted) bedGraph to disk ONCE instead
    # of write-unsorted → read → re-sort → rewrite. pipefail (set at top) makes
    # a genomecov failure fail the whole pipe.
    echo "TRACKS | COVERAGE | Computing positive strand (genomecov | sort)..."
    if ! bedtools genomecov \\
      -ibam "\${bam}" \\
      -\${end_type} \\
      -strand + \\
      -bg \\
      | LC_ALL=C sort -S "\${SORT_MEM}" -T "\${SORT_TMP}" --parallel="\${SORT_PARALLEL}" -k1,1 -k2,2n \\
      > "\${prefix}.pos.bedgraph"; then
      echo "TRACKS | ERROR | Failed to generate/sort positive strand coverage"
      return 1
    fi

    local pos_lines=\$(wc -l < "\${prefix}.pos.bedgraph" | tr -d ' ')
    local pos_size=\$(stat -c%s "\${prefix}.pos.bedgraph" 2>/dev/null || stat -f%z "\${prefix}.pos.bedgraph" 2>/dev/null || echo "unknown")
    echo "TRACKS | COVERAGE | Positive strand: \${pos_lines} regions (\${pos_size} bytes)"

    # Negative strand coverage (mirrored with -scale -1), sorted in one pass.
    echo "TRACKS | COVERAGE | Computing negative strand (mirrored with -scale -1 | sort)..."
    if ! bedtools genomecov \\
      -ibam "\${bam}" \\
      -\${end_type} \\
      -strand - \\
      -bg \\
      -scale -1 \\
      | LC_ALL=C sort -S "\${SORT_MEM}" -T "\${SORT_TMP}" --parallel="\${SORT_PARALLEL}" -k1,1 -k2,2n \\
      > "\${prefix}.neg.bedgraph"; then
      echo "TRACKS | ERROR | Failed to generate/sort negative strand coverage"
      return 1
    fi
    
    local neg_lines=\$(wc -l < "\${prefix}.neg.bedgraph" | tr -d ' ')
    local neg_size=\$(stat -c%s "\${prefix}.neg.bedgraph" 2>/dev/null || stat -f%z "\${prefix}.neg.bedgraph" 2>/dev/null || echo "unknown")
    echo "TRACKS | COVERAGE | Negative strand: \${neg_lines} regions (\${neg_size} bytes)"
    
    # Convert to BigWig
    echo "TRACKS | COVERAGE | Converting to BigWig format..."
    if ! make_bigwig "\${prefix}.pos.bedgraph" "\${prefix}.pos.bw"; then
      echo "TRACKS | ERROR | Failed to create positive strand BigWig"
      return 1
    fi
    if ! make_bigwig "\${prefix}.neg.bedgraph" "\${prefix}.neg.bw"; then
      echo "TRACKS | ERROR | Failed to create negative strand BigWig"
      return 1
    fi
    
    echo "TRACKS | COVERAGE | Complete: \${prefix}"
    return 0
  }

  # Perform UMI deduplication if enabled
  perform_umi_dedup() {
    local input_bam="\$1"
    local output_bam="\$2"
    local is_pe="\$3"
    local stats_file="\${4:-\${SAMPLE_ID}.dedup_stats.txt}"

    echo "TRACKS | DEDUP | Starting UMI deduplication..."
    echo "TRACKS | DEDUP | Mode: \$([ "\${is_pe}" == "true" ] && echo "Paired-end" || echo "Single-end")"
    
    # Parallelize across chromosomes: umi_tools calls duplicates from
    # (chromosome, position, UMI[, mate position]) only — it never compares
    # reads across chromosomes — so per-chromosome dedup is provably
    # identical to whole-BAM dedup. umi_tools has no thread flag of its own
    # and was observed pinned at ~95-98% of ONE core for 8-19min on
    # 91M/155M-read BAMs while the task holds \${THREADS} cpus (3 of 4 idle).
    # PE mates only land in the same chromosome shard when properly paired
    # (flag 2); everything else (discordant / mate-unmapped) dedups as ONE
    # unsplit job so cross-chromosome pairing is never broken. Shards are
    # extracted with samtools view from a coordinate-sorted BAM, so each
    # stays internally sorted and samtools merge k-way merges them straight
    # back into the original order.
    mkdir -p dedup_shards
    rm -f dedup_shards/* dedup_fail.flag 2>/dev/null || true

    mapfile -t DEDUP_CHROMS < <(samtools view -H "\${input_bam}" | \\
      awk -F'\\t' '\$1=="@SQ"{for(i=2;i<=NF;i++) if(\$i ~ /^SN:/){sub(/^SN:/,"",\$i); print \$i}}')

    dedup_shard() {
      local shard_in="\$1" shard_out="\$2" paired_flag="\$3"
      if ! umi_tools dedup \${paired_flag} -I "\${shard_in}" -S "\${shard_out}" \\
           --log="\${shard_out}.log" >>"\${shard_out}.err" 2>&1; then
        echo "FAIL: \${shard_in}" >> dedup_fail.flag
      fi
    }

    DEDUP_PIDS=()
    launch_dedup_shard() {
      while :; do
        local alive=0 p
        for p in "\${DEDUP_PIDS[@]:-}"; do
          [ -n "\${p}" ] && kill -0 "\${p}" 2>/dev/null && alive=\$(( alive + 1 ))
        done
        [ "\${alive}" -lt "\${THREADS}" ] && break
        sleep 0.5
      done
      dedup_shard "\$1" "\$2" "\$3" &
      DEDUP_PIDS+=(\$!)
    }

    SHARD_OUT_BAMS=()
    if [[ "\${is_pe}" == "true" ]]; then
      for chr in "\${DEDUP_CHROMS[@]}"; do
        safe_chr="\${chr//[^A-Za-z0-9_.-]/_}"
        shard_in="dedup_shards/pp_\${safe_chr}.bam"
        samtools view -@ "\${THREADS}" -b -f 2 -F 260 "\${input_bam}" "\${chr}" > "\${shard_in}"
        if [[ "\$(samtools view -c "\${shard_in}")" -eq 0 ]]; then rm -f "\${shard_in}"; continue; fi
        samtools index -@ "\${THREADS}" "\${shard_in}"
        shard_out="dedup_shards/pp_\${safe_chr}.dedup.bam"
        launch_dedup_shard "\${shard_in}" "\${shard_out}" "--paired"
        SHARD_OUT_BAMS+=("\${shard_out}")
      done
      other_in="dedup_shards/other.bam"
      samtools view -@ "\${THREADS}" -b -F 262 "\${input_bam}" > "\${other_in}"
      if [[ "\$(samtools view -c "\${other_in}")" -gt 0 ]]; then
        samtools index -@ "\${THREADS}" "\${other_in}"
        other_out="dedup_shards/other.dedup.bam"
        launch_dedup_shard "\${other_in}" "\${other_out}" "--paired"
        SHARD_OUT_BAMS+=("\${other_out}")
      fi
    else
      for chr in "\${DEDUP_CHROMS[@]}"; do
        safe_chr="\${chr//[^A-Za-z0-9_.-]/_}"
        shard_in="dedup_shards/\${safe_chr}.bam"
        samtools view -@ "\${THREADS}" -b "\${input_bam}" "\${chr}" > "\${shard_in}"
        if [[ "\$(samtools view -c "\${shard_in}")" -eq 0 ]]; then rm -f "\${shard_in}"; continue; fi
        samtools index -@ "\${THREADS}" "\${shard_in}"
        shard_out="dedup_shards/\${safe_chr}.dedup.bam"
        launch_dedup_shard "\${shard_in}" "\${shard_out}" ""
        SHARD_OUT_BAMS+=("\${shard_out}")
      done
    fi

    for p in "\${DEDUP_PIDS[@]}"; do
      wait "\${p}" || true
    done

    if [[ -s dedup_fail.flag ]]; then
      echo "TRACKS | ERROR | Failed dedup shard(s):"
      sed 's/^/TRACKS | ERROR |   /' dedup_fail.flag
      cat dedup_shards/*.err >&2 2>/dev/null || true
      echo "TRACKS | ERROR | umi_tools deduplication failed"
      return 1
    fi

    echo "TRACKS | DEDUP | Merging \${#SHARD_OUT_BAMS[@]} chromosome shard(s)..."
    if ! samtools merge -f -@ "\${THREADS}" "\${output_bam}" "\${SHARD_OUT_BAMS[@]}"; then
      echo "TRACKS | ERROR | Failed to merge deduplicated shards"
      return 1
    fi
    cat dedup_shards/*.dedup.bam.log > "\${stats_file}" 2>/dev/null || true
    rm -rf dedup_shards
    
    # Index deduplicated BAM
    with_io_lock samtools index -@ \${THREADS} "\${output_bam}"
    
    # Report statistics
    local before_reads=\$(samtools view -c -F \${FLAG_PRIMARY_MAPPED} "\${input_bam}")
    local after_reads=\$(samtools view -c -F \${FLAG_PRIMARY_MAPPED} "\${output_bam}")
    local removed=\$((before_reads - after_reads))
    local pct_removed=\$(awk -v b="\${before_reads}" -v r="\${removed}" 'BEGIN{printf "%.2f", (b>0)?(r*100.0/b):0}')
    
    echo "TRACKS | DEDUP | Reads before: \${before_reads}"
    echo "TRACKS | DEDUP | Reads after: \${after_reads}"
    echo "TRACKS | DEDUP | Removed: \${removed} (\${pct_removed}%)"
    
    # Append summary to stats file
    {
      echo ""
      echo "=== Summary ==="
      echo "reads_before=\${before_reads}"
      echo "reads_after=\${after_reads}"
      echo "reads_removed=\${removed}"
      echo "percent_removed=\${pct_removed}"
    } >> "\${stats_file}"
    
    return 0
  }

  ###########################################################################
  # 4) OPTIONAL UMI DEDUPLICATION
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | DEDUP | Checking UMI deduplication settings..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Start with main BAM
  INPUT_BAM="\${MAIN_BAM}"

  if [[ "\${UMI_ENABLED}" == "true" && \${UMI_LENGTH} -gt 0 ]]; then
    if command -v umi_tools >/dev/null 2>&1; then
      # Create working copy and index
      with_io_lock cp "\${MAIN_BAM}" aligned.bam
      with_io_lock samtools index -@ \${THREADS} aligned.bam
      
      # Perform deduplication
      if perform_umi_dedup "aligned.bam" "deduplicated.bam" "\${IS_PE}"; then
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
        } > "\${SAMPLE_ID}.dedup_stats.txt"
      fi
    else
      echo "TRACKS | WARNING | umi_tools not found, skipping deduplication"
      {
        echo "umi_deduplication=skipped"
        echo "status=umi_tools not available"
        echo "reads_before=N/A"
        echo "reads_after=N/A"
        echo "reads_removed=N/A"
      } > "\${SAMPLE_ID}.dedup_stats.txt"
    fi
  else
    echo "TRACKS | DEDUP | UMI deduplication disabled"
    {
      echo "umi_deduplication=disabled"
      echo "status=not requested in parameters"
      echo "reads_before=N/A"
      echo "reads_after=N/A"
      echo "reads_removed=N/A"
    } > "\${SAMPLE_ID}.dedup_stats.txt"
  fi

  # ── Optional: also UMI-deduplicate the allMap BAM (off by default) ──────────
  # When UMI dedup is on, the MAIN tracks use the deduped BAM but allMap tracks
  # use the raw allMap BAM, so the two are deduped inconsistently. Enable
  # params.umi.dedup_allmap=true to dedup allMap too. OFF by default because
  # umi_tools dedup is designed for unique alignments; on a bowtie2 -k multimapper
  # BAM (NH>1, MAPQ=255) its position-based dedup is approximate. Falls back to the
  # original allMap BAM on any failure, so the standard path is never broken.
  DEDUP_ALLMAP="${(params.umi?.dedup_allmap == true) ? 'true' : 'false'}"
  if [[ "\${UMI_ENABLED}" == "true" && \${UMI_LENGTH} -gt 0 && "\${DEDUP_ALLMAP}" == "true" ]] \\
     && command -v umi_tools >/dev/null 2>&1; then
    echo "TRACKS | DEDUP | Also UMI-deduplicating allMap BAM (experimental for multimappers)..."
    with_io_lock cp "\${ALLMAP_BAM}" allmap_in.bam
    with_io_lock samtools index -@ \${THREADS} allmap_in.bam
    if perform_umi_dedup "allmap_in.bam" "allmap_dedup.bam" "\${IS_PE}" "allmap_dedup_stats.txt"; then
      ALLMAP_BAM="allmap_dedup.bam"
      echo "TRACKS | DEDUP | allMap deduplicated"
    else
      echo "TRACKS | WARNING | allMap dedup failed; using original allMap BAM"
    fi
  fi

  # NOTE: bam_for_downstream.bam (consumed by module 11 for Pol-II gene metrics)
  # is created AFTER the PE mate-filtering step below, so that in paired-end mode
  # it is the SAME signal-mate-only BAM the coverage tracks use. Previously it was
  # the full (both-mate) deduped BAM, so gene metrics counted the noise mate too,
  # inflating tss_cpm/body_cpm relative to the tracks. See section 4c.

  ###########################################################################
  # 4b) PE MATE FILTERING FOR 3'/5' COVERAGE TRACKS
  ###########################################################################
  #
  # PRO-seq PE alignment layout (from module 05):
  #   bowtie2 --ff  -1 R2  -2 RC(R1)
  # In the BAM this becomes:
  #   Read1 (flag 64)  = original R2  → 5' end of fragment, NOT the Pol II position
  #   Read2 (flag 128) = RC(R1)       → 3' end of nascent RNA = Pol II position ✓
  #
  # Using the full paired BAM for -3 coverage would add one noise hit (from the
  # non-signal mate) for every correct Pol II hit, distorting track shapes. For
  # SE data both reads carry signal, so no filtering is needed.
  #
  # Which mate carries the signal end is configurable: params.align.pe_signal_mate
  # = 'read2' (default, flag 128 = RC(R1), matches the module-05 layout above) or
  # 'read1' (flag 64) for chemistries where R1 carries the signal end. NOTE: using
  # 'read1' generally also requires a matching alignment orientation in module 05.
  #
  # The allMap BAM is mate-filtered the same way so allMap tracks stay comparable.
  PE_SIGNAL_MATE="${params.align?.pe_signal_mate ?: 'read2'}"
  if [[ "\${PE_SIGNAL_MATE}" == "read1" ]]; then PE_MATE_FLAG=64; else PE_MATE_FLAG=128; fi

  BAM_FOR_COVERAGE="\${INPUT_BAM}"
  ALLMAP_BAM_FOR_COVERAGE="\${ALLMAP_BAM}"

  if [[ "\${IS_PE}" == "true" ]]; then
    echo "────────────────────────────────────────────────────────────────────────"
    echo "TRACKS | PE_FILTER | Paired-end: keeping only the signal mate (\${PE_SIGNAL_MATE}, flag \${PE_MATE_FLAG}) for coverage..."
    echo "────────────────────────────────────────────────────────────────────────"

    samtools view -@ "\${THREADS}" -f \${PE_MATE_FLAG} -b "\${INPUT_BAM}" \\
      | samtools sort -@ "\${THREADS}" -o pe_r2_main.bam
    samtools index -@ "\${THREADS}" pe_r2_main.bam

    samtools view -@ "\${THREADS}" -f \${PE_MATE_FLAG} -b "\${ALLMAP_BAM}" \\
      | samtools sort -@ "\${THREADS}" -o pe_r2_allmap.bam
    samtools index -@ "\${THREADS}" pe_r2_allmap.bam

    R2_MAIN_COUNT=\$(samtools view -c -F 4 pe_r2_main.bam)
    R2_ALLMAP_COUNT=\$(samtools view -c -F 4 pe_r2_allmap.bam)
    echo "TRACKS | PE_FILTER | Main BAM signal-mate count:   \${R2_MAIN_COUNT}"
    echo "TRACKS | PE_FILTER | AllMap BAM signal-mate count: \${R2_ALLMAP_COUNT}"

    BAM_FOR_COVERAGE="pe_r2_main.bam"
    ALLMAP_BAM_FOR_COVERAGE="pe_r2_allmap.bam"
    echo "TRACKS | PE_FILTER | Coverage tracks will use signal-mate-only BAMs"
  else
    echo "TRACKS | PE_FILTER | Single-end mode: using full BAM for coverage"
  fi

  ###########################################################################
  # 4c) BAM HANDED TO DOWNSTREAM POL-II METRICS (module 11)
  ###########################################################################
  #
  # Use the SAME BAM the coverage tracks use: in PE this is the Read2-only
  # (flag 128 = RC(R1) = Pol II 3' end) BAM, in SE the full deduped BAM. This
  # keeps gene-level TSS/body counts consistent with the published tracks
  # instead of counting both mates in PE.
  echo "TRACKS | OUTPUT | Copying BAM used for tracks + Pol-II metrics (deduped when UMI on)..."
  with_io_lock cp "\${BAM_FOR_COVERAGE}" bam_for_downstream.bam
  with_io_lock samtools index -@ \${THREADS} bam_for_downstream.bam

  ###########################################################################
  # 5) PREPARE GENOME SIZES
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | GENOME | Extracting chromosome sizes..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Get chromosome sizes from BAM header (more reliable than FASTA).
  # Use BAM_FOR_COVERAGE so the sizes match the BAM we'll actually feed to bedtools.
  samtools view -H "\${BAM_FOR_COVERAGE}" | \\
    grep '^@SQ' | \\
    cut -f2,3 | \\
    sed 's/SN://g' | \\
    sed 's/LN://g' \\
    > genome.sizes

  if [[ ! -s genome.sizes ]]; then
    tracktx_error "generate_coverage_tracks" "Failed to extract chromosome sizes from BAM header" "Check BAM file integrity"
  fi

  CHR_COUNT=\$(wc -l < genome.sizes | tr -d ' ')
  TOTAL_SIZE=\$(awk '{sum+=\$2} END{print sum}' genome.sizes)
  echo "TRACKS | GENOME | Chromosomes: \${CHR_COUNT}"
  echo "TRACKS | GENOME | Total genome size: \${TOTAL_SIZE} bp"

  # Display first few chromosomes
  echo "TRACKS | GENOME | First chromosomes:"
  head -5 genome.sizes | sed 's/^/TRACKS | GENOME |   /'

  ###########################################################################
  # 6) GENERATE 3' END COVERAGE (Always)
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | 3P | Generating 3' end coverage tracks..."
  echo "────────────────────────────────────────────────────────────────────────"

  # The four coverage jobs are independent, but each is NOT free in RAM:
  #   • bedtools genomecov -ibam allocates a full-length counts array for the
  #     current chromosome (~1 GB for hs1 chr1), and
  #   • its sort buffers up to SORT_MEM.
  # Running all four at once (genomecov + sort × 4) overruns the task memory
  # cgroup on memory-modest hosts and the kernel OOM-kills a genomecov
  # (observed: allMap 3' "Killed" with this process capped at a few GB). So cap
  # concurrency by the task's memory budget instead of always launching 4.
  #
  # Budget each running job at ~2 GB peak (genomecov array + sort buffer) and
  # use up to ~70% of the task memory, never more than the allotted CPUs.
  MEM_GB=${task.memory.toGiga()}
  MAX_PAR=\$(( MEM_GB * 70 / 100 / 2 ))
  [ "\${MAX_PAR}" -lt 1 ] && MAX_PAR=1
  [ "\${MAX_PAR}" -gt ${task.cpus} ] && MAX_PAR=${task.cpus}
  MAX_PAR=\${TRACKS_MAX_PAR:-\${MAX_PAR}}

  # Per-sort memory scaled to the chosen concurrency (≤50% of budget shared
  # across the at-most-MAX_PAR concurrent sorts). Temp on a fast dir, not USB.
  SORT_MEM=\$(( MEM_GB * 50 / 100 / MAX_PAR ))
  [ "\${SORT_MEM}" -lt 1 ] && SORT_MEM=1
  export SORT_MEM="\${SORT_MEM}G"
  export SORT_TMP="\${SORT_TMPDIR:-/tmp}"
  # Docker gets --cpus=THREADS (a CFS quota) but NOT --cpuset-cpus, so
  # /proc/cpuinfo inside the container still reports the HOST's full core
  # count. GNU sort's --parallel auto-detection reads that and may spawn
  # threads for all host cores while only actually getting THREADS worth of
  # scheduling — with MAX_PAR of these sorts running at once, that's
  # threads-per-sort × MAX_PAR competing for THREADS cpu-shares, i.e.
  # self-inflicted contention, not speed. Cap each sort's own parallelism at
  # its fair share of THREADS instead of letting it guess wrong.
  SORT_PARALLEL=\$(( THREADS / MAX_PAR ))
  [ "\${SORT_PARALLEL}" -lt 1 ] && SORT_PARALLEL=1
  export SORT_PARALLEL
  echo "TRACKS | 3P | Concurrency: \${MAX_PAR} parallel job(s), SORT_MEM=\${SORT_MEM} each, SORT_PARALLEL=\${SORT_PARALLEL} (task mem=\${MEM_GB}G)"

  # Throttled fan-out: never let more than MAX_PAR generate_coverage run at once.
  # Each job records its own failure to a flag file so all per-job errors land in
  # the log before we abort. Log output is interleaved; each job identifies
  # itself via its BAM filename in the TRACKS | COVERAGE messages.
  #
  # IMPORTANT: track our OWN job PIDs explicitly and wait on those PIDs only.
  # We must NOT use bare \`wait\` or \`wait -n\`: the \`exec > >(tee ...)\` process
  # substitutions at the top of this script spawn long-lived \`tee\` children that
  # never exit until the script's FDs close, so a bare wait would block on them
  # forever ("process hasn't exited" → Nextflow kills the hung task). Likewise we
  # throttle by polling our tracked PIDs with \`kill -0\`, not \`jobs\`.
  COV_FAIL_FLAG="cov_fail.flag"
  rm -f "\${COV_FAIL_FLAG}"

  run_cov() {
    generate_coverage "\$1" "\$2" "\$3" || echo "FAIL: \$3" >> "\${COV_FAIL_FLAG}"
  }

  COV_PIDS=()
  launch_cov() {
    # Block until fewer than MAX_PAR of OUR jobs are still alive.
    while :; do
      local alive=0 p
      for p in "\${COV_PIDS[@]:-}"; do
        [ -n "\${p}" ] && kill -0 "\${p}" 2>/dev/null && alive=\$(( alive + 1 ))
      done
      [ "\${alive}" -lt "\${MAX_PAR}" ] && break
      sleep 0.5
    done
    run_cov "\$1" "\$2" "\$3" &
    COV_PIDS+=(\$!)
  }

  # Use the PE-filtered BAMs (Read2-only) in PE mode; full BAMs in SE mode.
  launch_cov "\${BAM_FOR_COVERAGE}"        "3" "3p/\${SAMPLE_ID}.3p"
  launch_cov "\${ALLMAP_BAM_FOR_COVERAGE}" "3" "3p/\${SAMPLE_ID}.allMap.3p"
  launch_cov "\${BAM_FOR_COVERAGE}"        "5" "5p/\${SAMPLE_ID}.5p"
  launch_cov "\${ALLMAP_BAM_FOR_COVERAGE}" "5" "5p/\${SAMPLE_ID}.allMap.5p"

  # Wait on each of our jobs specifically (never bare wait — see note above).
  for p in "\${COV_PIDS[@]}"; do
    wait "\${p}" || true
  done

  if [[ -s "\${COV_FAIL_FLAG}" ]]; then
    echo "TRACKS | ERROR | Failed coverage jobs:"
    sed 's/^/TRACKS | ERROR |   /' "\${COV_FAIL_FLAG}"
    tracktx_error "generate_coverage_tracks" "One or more coverage generation jobs failed" "Check tracks.log for per-job error messages"
  fi

  # NOTE: All four track sets (main/allMap × 3'/5') are produced by the single
  # parallel block above, using the PE mate-filtered BAMs (Read2-only) in PE
  # mode and the full BAMs in SE mode. The previous sequential re-runs here
  # regenerated allMap-3', main-5' and allMap-5' from the *unfiltered* BAMs,
  # overwriting the correct Read2-only outputs in PE mode (contaminating tracks
  # with the wrong mate's end position). They were removed so the mate-filtered
  # outputs survive for paired-end data.

  echo "TRACKS | 3P | 3' end coverage complete"
  echo "TRACKS | 5P | 5' end coverage complete"

  ###########################################################################
  # 8) CREATE DOCUMENTATION
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | README | Creating documentation..."
  echo "────────────────────────────────────────────────────────────────────────"

  cat > \${SAMPLE_ID}.README_tracks.txt <<DOCEOF
================================================================================
COVERAGE TRACKS — ${sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Strand-specific coverage tracks at nucleotide resolution for PRO-seq analysis.
  
  Track Types:
    • 3' end coverage: Always generated (PRO-seq standard)
    • 5' end coverage: Always generated (PE and SE)
  
  BAM Sources:
    • Main BAM: Primary alignments ${params.umi?.enabled ? 'with UMI deduplication' : '(duplicates retained)'}
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
  3p/${sample_id}.3p.pos.bedgraph    — Positive strand (main BAM)
  3p/${sample_id}.3p.neg.bedgraph    — Negative strand (main BAM, mirrored)
  3p/${sample_id}.3p.pos.bw          — BigWig format (positive)
  3p/${sample_id}.3p.neg.bw          — BigWig format (negative)
  
  3p/${sample_id}.allMap.3p.pos.bedgraph — AllMap BAM positive strand
  3p/${sample_id}.allMap.3p.neg.bedgraph — AllMap BAM negative strand (mirrored)
  3p/${sample_id}.allMap.3p.pos.bw       — BigWig format
  3p/${sample_id}.allMap.3p.neg.bw       — BigWig format

5' End Coverage:
  5p/${sample_id}.5p.pos.bedgraph        — Main BAM positive strand
  5p/${sample_id}.5p.neg.bedgraph        — Main BAM negative strand (mirrored)
  5p/${sample_id}.5p.pos.bw              — BigWig format
  5p/${sample_id}.5p.neg.bw              — BigWig format
  
  5p/${sample_id}.allMap.5p.pos.bedgraph — AllMap BAM positive strand
  5p/${sample_id}.allMap.5p.neg.bedgraph — AllMap BAM negative strand (mirrored)
  5p/${sample_id}.allMap.5p.pos.bw       — BigWig format
  5p/${sample_id}.allMap.5p.neg.bw       — BigWig format

Statistics:
  ${sample_id}.dedup_stats.txt       — UMI deduplication statistics
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
  Status: ${params.umi?.enabled ? 'Enabled' : 'Disabled'}
  ${params.umi?.enabled ? 'Length: ' + params.umi.length + ' bp' : ''}
  ${params.umi?.enabled ? 'Duplicates removed before track generation' : 'Duplicates retained in coverage'}

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
    • NOTE: divergent transcription detection uses the MAIN (primary) 3' tracks,
      not allMap (see main.nf STEP 10 wiring)

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These raw tracks will be:
  1. Normalized to CPM and siCPM (next module)
  2. Used for divergent transcription detection (main 3' tracks)
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
  • 16 bedGraph + 16 BigWig files (32 total) for both PE and SE
    - 8 files for 3' end coverage (4 main + 4 allMap)
    - 8 files for 5' end coverage (4 main + 4 allMap)

Troubleshooting:
  • Empty bedGraph: Check if BAM has mapped reads
  • BigWig conversion failure: bedGraph sorting issue (automatically handled)
  • Large file sizes: Expected for high-coverage samples

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • All coverage values are raw counts (not normalized)
  • Negative strand uses -scale -1 (native bedtools feature)
  • Chromosome sizes from BAM header ensure coordinate consistency
  • BigWig creation uses 600-second timeout for large genomes
  • bedGraph sorting is mandatory (automatically performed)

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  UMI deduplication:    ${params.umi?.enabled ? 'Enabled' : 'Disabled'}
  UMI length:           ${params.umi?.length ?: 'N/A'} bp
  Library type:         ${is_paired == "true" ? "Paired-end" : "Single-end"}
  CPU threads:          ${task.cpus}

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Module:   06_generate_coverage_tracks
  Date:     \$(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample:   ${sample_id}

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
  BG_COUNT=\$(find 3p 5p -name "*.bedgraph" -type f 2>/dev/null | wc -l | tr -d ' ')
  BW_COUNT=\$(find 3p 5p -name "*.bw" -type f 2>/dev/null | wc -l | tr -d ' ')

  # Expected file counts (8 for 3p + 8 for 5p, both PE and SE)
  EXPECTED_BG=16
  EXPECTED_BW=16

  echo "TRACKS | VALIDATE | bedGraph files: \${BG_COUNT}/\${EXPECTED_BG}"
  echo "TRACKS | VALIDATE | BigWig files: \${BW_COUNT}/\${EXPECTED_BW}"

  # Check critical 3' files (must be non-empty)
  for file in \\
    "3p/\${SAMPLE_ID}.3p.pos.bedgraph" \\
    "3p/\${SAMPLE_ID}.3p.neg.bedgraph" \\
    "3p/\${SAMPLE_ID}.3p.pos.bw" \\
    "3p/\${SAMPLE_ID}.3p.neg.bw" \\
    "3p/\${SAMPLE_ID}.allMap.3p.pos.bedgraph" \\
    "3p/\${SAMPLE_ID}.allMap.3p.neg.bedgraph" \\
    "3p/\${SAMPLE_ID}.allMap.3p.pos.bw" \\
    "3p/\${SAMPLE_ID}.allMap.3p.neg.bw"; do
    
    if [[ ! -s "\${file}" ]]; then
      tracktx_error "generate_coverage_tracks" "Missing or empty critical file: \${file}" "Check tracks.log in work dir"
    else
      FILE_SIZE=\$(stat -c%s "\${file}" 2>/dev/null || stat -f%z "\${file}" 2>/dev/null || echo "unknown")
      echo "TRACKS | VALIDATE | \${file}: \${FILE_SIZE} bytes"
    fi
  done

  # Check that 5' files exist
  for file in \\
    "5p/\${SAMPLE_ID}.5p.pos.bedgraph" \\
    "5p/\${SAMPLE_ID}.5p.neg.bedgraph" \\
    "5p/\${SAMPLE_ID}.5p.pos.bw" \\
    "5p/\${SAMPLE_ID}.5p.neg.bw" \\
    "5p/\${SAMPLE_ID}.allMap.5p.pos.bedgraph" \\
    "5p/\${SAMPLE_ID}.allMap.5p.neg.bedgraph" \\
    "5p/\${SAMPLE_ID}.allMap.5p.pos.bw" \\
    "5p/\${SAMPLE_ID}.allMap.5p.neg.bw"; do
    
    if [[ ! -f "\${file}" ]]; then
      tracktx_error "generate_coverage_tracks" "Missing file: \${file}" "Check tracks.log in work dir"
    fi
  done


  echo "TRACKS | VALIDATE | All required files present"

  ###########################################################################
  # FINAL SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | SUMMARY | Processing Complete"
  echo "────────────────────────────────────────────────────────────────────────"
  echo "TRACKS | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "TRACKS | SUMMARY | Library type: \$([ "\${IS_PE}" == "true" ] && echo "Paired-end" || echo "Single-end") (3' + 5')"
  echo "TRACKS | SUMMARY | UMI dedup: \${UMI_ENABLED}"
  echo "TRACKS | SUMMARY | bedGraph files: \${BG_COUNT}"
  echo "TRACKS | SUMMARY | BigWig files: \${BW_COUNT}"
  echo "TRACKS | SUMMARY | Total output size: \$(du -sh . 2>/dev/null | cut -f1 || echo "unknown")"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "TRACKS | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}
