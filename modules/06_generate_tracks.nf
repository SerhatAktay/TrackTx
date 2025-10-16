// ============================================================================
// generate_tracks.nf — 3′/5′ stranded coverage (+ allMap), per-sample
// ----------------------------------------------------------------------------
// Overview
//   • (Opt) UMI dedup → BAM (+ .bai) used for coverage
//   • BAM → BED → bedGraph pipeline with genome bounds validation (-g flag)
//   • Main BAM: 3′ POS/NEG (NEG mirrored); PE also 5′ POS/NEG
//   • allMap BAM: same (3′ always; 5′ if PE)
//   • BigWigs via UCSC bedGraphToBigWig
//   • Unified step log; README
//
// Why it's cleaner
//   • bamToBed | bedtools genomecov -g ensures valid coordinates upfront
//   • No post-hoc clipping needed (genome.sizes prevents out-of-bounds)
//   • Streams through pipes (no large intermediate files written to disk)
//
// Inputs
//   tuple( sample_id, filtered_bam, spikein_bam, condition, timepoint, replicate )
//   genome_fa (FASTA; .fai will be created if missing)
//   is_paired  ("true" | "false")
//   allmap_bam (sample_allMap.bam; coord-sorted)
//
// Outputs (publishDir):
//   ${params.output_dir}/03_genome_tracks/${sample_id}/
//     ├── 3p/${sample_id}.3p.pos.bedgraph / .bw
//     ├── 3p/${sample_id}.3p.neg.bedgraph / .bw
//     ├── 3p/${sample_id}.allMap.3p.pos.bedgraph / .bw
//     ├── 3p/${sample_id}.allMap.3p.neg.bedgraph / .bw
//     ├── 5p/${sample_id}.5p.pos.bedgraph / .bw                 (PE only)
//     ├── 5p/${sample_id}.5p.neg.bedgraph / .bw                 (PE only)
//     ├── 5p/${sample_id}.allMap.5p.pos.bedgraph / .bw          (PE only)
//     ├── 5p/${sample_id}.allMap.5p.neg.bedgraph / .bw          (PE only)
//     ├── ${sample_id}.dedup_stats.txt
//     ├── ${sample_id}.README_tracks.txt
//     └── tracks.log
// ============================================================================

nextflow.enable.dsl = 2

process generate_tracks {

  // ── Meta / resources ─────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/03_genome_tracks/${sample_id}", mode: 'copy', overwrite: true

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id), path(filtered_bam), path(spikein_bam),
          val(condition), val(timepoint), val(replicate)
    path genome_fa
    val  is_paired
    path allmap_bam

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    // Main 3′ (always)
    tuple val(sample_id),
          path("3p/${sample_id}.3p.pos.bedgraph"),
          path("3p/${sample_id}.3p.neg.bedgraph"),
          path("3p/${sample_id}.3p.pos.bw"),
          path("3p/${sample_id}.3p.neg.bw"),
          val(condition), val(timepoint), val(replicate),
          emit: bw3p_pair

    // Main 5′ (PE only)
    tuple val(sample_id),
          path("5p/${sample_id}.5p.pos.bedgraph", optional: true),
          path("5p/${sample_id}.5p.neg.bedgraph", optional: true),
          path("5p/${sample_id}.5p.pos.bw",       optional: true),
          path("5p/${sample_id}.5p.neg.bw",       optional: true),
          val(condition), val(timepoint), val(replicate),
          emit: bw5p_pair

    // allMap 3′ (always)
    tuple val(sample_id),
          path("3p/${sample_id}.allMap.3p.pos.bedgraph"),
          path("3p/${sample_id}.allMap.3p.neg.bedgraph"),
          path("3p/${sample_id}.allMap.3p.pos.bw"),
          path("3p/${sample_id}.allMap.3p.neg.bw"),
          val(condition), val(timepoint), val(replicate),
          emit: allmap3p_pair

    // allMap 5′ (PE only)
    tuple val(sample_id),
          path("5p/${sample_id}.allMap.5p.pos.bedgraph", optional: true),
          path("5p/${sample_id}.allMap.5p.neg.bedgraph", optional: true),
          path("5p/${sample_id}.allMap.5p.pos.bw",       optional: true),
          path("5p/${sample_id}.allMap.5p.neg.bw",       optional: true),
          val(condition), val(timepoint), val(replicate),
          emit: allmap5p_pair

    // Legacy convenience tuple
    tuple val(sample_id),
          path(filtered_bam), path(spikein_bam),
          path("3p/${sample_id}.3p.pos.bedgraph"),
          path("3p/${sample_id}.3p.neg.bedgraph"),
          path("5p/${sample_id}.5p.pos.bedgraph", optional: true),
          path("5p/${sample_id}.5p.neg.bedgraph", optional: true),
          val(condition), val(timepoint), val(replicate),
          emit: track_tuple

    // Logs / docs
    tuple val(sample_id), path("${sample_id}.dedup_stats.txt"),
          val(condition), val(timepoint), val(replicate), emit: dedup_stats
    path "${sample_id}.README_tracks.txt"
    path "tracks.log", emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a tracks.log) 2>&1
  trap 'echo "ERROR  [tracks] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR

  # ── Bindings / locals ────────────────────────────────────────────────────
  THREADS=!{task.cpus}
  PFX="!{sample_id}"
  BAM="!{filtered_bam}"            # main sample BAM (coord-sorted)
  ALLMAP_BAM="!{allmap_bam}"       # all mapped BAM (coord-sorted)
  IS_PE="!{is_paired}"
  UMI_ENABLED="!{ params.umi?.enabled ? 'true' : 'false' }"
  UMI_LEN=!{ params.umi?.length ?: 0 }
  FORCE_SORT_BG="!{ params.force_sort_bedgraph ? 'true' : 'false' }"

  echo "INFO  [tracks] ▶ sample=${PFX}  PE=${IS_PE}  cpus=${THREADS} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  [[ -s "$BAM" ]] || { echo "ERROR  Empty filtered BAM" >&2; exit 1; }
  [[ -s "$ALLMAP_BAM" ]] || { echo "ERROR  Missing allMap BAM — required for 8× 3′ outputs" >&2; exit 1; }

  command -v samtools         >/dev/null || { echo "ERROR samtools not found" >&2; exit 1; }
  command -v bedtools         >/dev/null || { echo "ERROR bedtools not found" >&2; exit 1; }
  command -v bedGraphToBigWig >/dev/null || { echo "ERROR bedGraphToBigWig not found" >&2; exit 1; }
  # bedClip optional

  # ─────────────────────────────────────────────────────────────────────────
  # 1) Optional UMI deduplication
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  Step 1/6: UMI deduplication (if enabled)"
  ln -sf "$BAM" aligned.bam
  samtools index -@ ${THREADS} aligned.bam
  INPUT_BAM=aligned.bam

    if [[ "$UMI_ENABLED" == "true" && ${UMI_LEN} -gt 0 ]]; then
      if command -v umi_tools >/dev/null 2>&1; then
        echo "INFO  … running umi_tools dedup (len=${UMI_LEN})"
        if [[ "${IS_PE}" == "true" ]]; then
          umi_tools dedup --paired -I aligned.bam -S deduplicated.bam --log=${PFX}.dedup_stats.txt
        else
          umi_tools dedup -I aligned.bam -S deduplicated.bam --log=${PFX}.dedup_stats.txt
        fi
        INPUT_BAM=deduplicated.bam
        samtools index -@ ${THREADS} "$INPUT_BAM"
      else
        echo "WARN  umi_tools not found — skipping UMI dedup" | tee ${PFX}.dedup_stats.txt
      fi
    else
      echo "UMI disabled" > ${PFX}.dedup_stats.txt
    fi

  # ─────────────────────────────────────────────────────────────────────────
  # 1.5) Disk space check before intensive operations
  # ─────────────────────────────────────────────────────────────────────────
  AVAIL_GB=$(df -BG . 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G' || echo "0")
  if [[ "$AVAIL_GB" -lt 5 ]]; then
    echo "WARNING [tracks] Low disk space: ${AVAIL_GB}GB available" >&2
    echo "WARNING [tracks] Track generation may fail with large BAMs" >&2
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 2) Genome sizes (for BigWig)
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  Step 2/6: preparing genome.sizes"
  if [[ -s "!{genome_fa}.fai" ]]; then
    cut -f1,2 "!{genome_fa}.fai" > genome.sizes
  else
    samtools faidx "!{genome_fa}" >/dev/null 2>&1
    cut -f1,2 "!{genome_fa}.fai" > genome.sizes
  fi
  [[ -s genome.sizes ]] || { echo "ERROR  Failed to build genome.sizes" >&2; exit 1; }

  # ─────────────────────────────────────────────────────────────────────────
  # 3) Helpers
  # ─────────────────────────────────────────────────────────────────────────
  make_bw () {   # bedGraph → BigWig (optional sort)
    local bg="$1" bw="$2"
    [[ -s "$bg" ]] || return 0
    
    echo "INFO  [tracks] Converting ${bg} to BigWig format"
    if [[ "$FORCE_SORT_BG" == "true" ]]; then
      echo "INFO  [tracks] Sorting ${bg} before BigWig conversion"
      LC_ALL=C sort -k1,1 -k2,2n -o "$bg" "$bg"
    fi
    
    echo "INFO  [tracks] Running bedGraphToBigWig on ${bg}"
    echo "INFO  [tracks] Input bedGraph size: $(wc -l < "$bg") lines"
    echo "INFO  [tracks] Genome sizes file: $(cat genome.sizes | wc -l) chromosomes"
    
    if ! timeout 600 bedGraphToBigWig "$bg" genome.sizes "$bw"; then
      echo "ERROR  [tracks] bedGraphToBigWig failed or timed out for $bg" >&2
      echo "ERROR  [tracks] Check if bedGraph file is valid and genome.sizes is correct" >&2
      return 1
    fi
    echo "INFO  [tracks] Successfully created ${bw} ($(stat -f%z "$bw" 2>/dev/null || stat -c%s "$bw" 2>/dev/null || echo "unknown") bytes)"
  }

  cov_pair () {  # stranded coverage for end type (3 or 5); NEG mirrored
    local bam="$1" endflag="$2" prefix="$3"
    echo "INFO  [tracks] Processing BAM: $bam for ${endflag}' end"
    
    # Check if BAM is valid before processing
    if ! samtools quickcheck "$bam" 2>/dev/null; then
      echo "ERROR  [tracks] Invalid BAM file: $bam" >&2
      return 1
    fi
    
    # POS strand: BAM → BED → bedGraph (with genome bounds validation)
    echo "INFO  [tracks] Computing + strand coverage for ${endflag}' end (via BED intermediate)"
    bamToBed -i "$bam" | bedtools genomecov -bg -i stdin -g genome.sizes -${endflag} -strand + > "${prefix}.pos.bedgraph"
    if [ $? -ne 0 ]; then
      echo "ERROR  [tracks] Failed to compute + strand coverage" >&2
      return 1
    fi
    
    # NEG strand (mirrored): BAM → BED → bedGraph (with genome bounds validation)
    echo "INFO  [tracks] Computing - strand coverage for ${endflag}' end (via BED intermediate, mirrored)"
    bamToBed -i "$bam" | bedtools genomecov -bg -i stdin -g genome.sizes -${endflag} -strand - | awk 'BEGIN{OFS="\\t"}{\$4=-\$4; print}' > "${prefix}.neg.bedgraph"
    if [ $? -ne 0 ]; then
      echo "ERROR  [tracks] Failed to compute - strand coverage" >&2
      return 1
    fi
    
    # No clipping needed - bedtools genomecov with -g flag handles bounds automatically!
    echo "INFO  [tracks] bedGraph files validated via -g genome.sizes (no clipping needed)"
    
    # Convert to BigWig
    echo "INFO  [tracks] Converting to BigWig format"
    make_bw "${prefix}.pos.bedgraph" "${prefix}.pos.bw"
    make_bw "${prefix}.neg.bedgraph" "${prefix}.neg.bw"
  }

  mkdir -p 3p 5p

  # ─────────────────────────────────────────────────────────────────────────
  # 4) Main 3′/5′ coverage (sequential to avoid resource contention)
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  Step 3/6: main BAM coverage → 3′"
  if ! cov_pair "$INPUT_BAM" 3 "3p/${PFX}.3p"; then
    echo "ERROR  [tracks] Failed to process main BAM 3′ coverage" >&2
    exit 1
  fi
  
  if [[ "${IS_PE}" == "true" ]]; then
    echo "INFO  Step 4/6: main BAM coverage → 5′ (PE)"
    if ! cov_pair "$INPUT_BAM" 5 "5p/${PFX}.5p"; then
      echo "ERROR  [tracks] Failed to process main BAM 5′ coverage" >&2
      exit 1
    fi
  fi
  echo "INFO  … main BAM coverage done"

  # ─────────────────────────────────────────────────────────────────────────
  # 5) allMap 3′/5′ coverage (sequential to avoid resource contention)
  # ─────────────────────────────────────────────────────────────────────────
  echo "INFO  Step 5/6: allMap BAM coverage → 3′"
  if ! cov_pair "$ALLMAP_BAM" 3 "3p/${PFX}.allMap.3p"; then
    echo "ERROR  [tracks] Failed to process allMap BAM 3′ coverage" >&2
    exit 1
  fi
  
  if [[ "${IS_PE}" == "true" ]]; then
    echo "INFO  Step 6/6: allMap BAM coverage → 5′ (PE)"
    if ! cov_pair "$ALLMAP_BAM" 5 "5p/${PFX}.allMap.5p"; then
      echo "ERROR  [tracks] Failed to process allMap BAM 5′ coverage" >&2
      exit 1
    fi
  fi
  echo "INFO  … allMap BAM coverage done"

  # ─────────────────────────────────────────────────────────────────────────
  # 5′ placeholders for SE libraries (to satisfy optional outputs)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${IS_PE}" != "true" ]]; then
    echo "INFO  [tracks] PE=false; creating placeholder 5′ outputs"
    : > "5p/${PFX}.5p.pos.bedgraph"
    : > "5p/${PFX}.5p.neg.bedgraph"
    : > "5p/${PFX}.5p.pos.bw"
    : > "5p/${PFX}.5p.neg.bw"
    : > "5p/${PFX}.allMap.5p.pos.bedgraph"
    : > "5p/${PFX}.allMap.5p.neg.bedgraph"
    : > "5p/${PFX}.allMap.5p.pos.bw"
    : > "5p/${PFX}.allMap.5p.neg.bw"
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 6) README
  # ─────────────────────────────────────────────────────────────────────────
  cat > ${PFX}.README_tracks.txt <<TXT
Genome tracks for sample: !{sample_id}

Main (from INPUT_BAM${UMI_ENABLED=='true' && UMI_LEN>0 ? ' after UMI dedup' : ''}):
  3′: 3p/${PFX}.3p.pos/neg.(bedgraph|bw)  (NEG mirrored)
  5′: 5p/${PFX}.5p.pos/neg.(bedgraph|bw)  (PE only)

allMap (from sample_allMap.bam):
  3′: 3p/${PFX}.allMap.3p.pos/neg.(bedgraph|bw)
  5′: 5p/${PFX}.allMap.5p.pos/neg.(bedgraph|bw)  (PE only)

Processing pipeline
  1. BAM → BED conversion (bamToBed)
  2. BED → bedGraph with genome bounds validation (bedtools genomecov -g)
  3. bedGraph → BigWig (bedGraphToBigWig)

Notes
  - bamToBed streams into bedtools genomecov for efficient processing
  - -g genome.sizes flag ensures all coordinates are valid (no clipping needed)
  - -3 / -5 flags extract end-specific coverage (single-nucleotide resolution)
  - Negative-strand tracks are mirrored (score *= -1) for UCSC-style viewing
  - BigWigs built with bedGraphToBigWig; set params.force_sort_bedgraph=true
    if UCSC sorting complaints appear on your system
TXT

  echo "INFO  [tracks] ✔ done for ${PFX} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
