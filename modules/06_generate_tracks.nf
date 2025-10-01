// ============================================================================
// generate_tracks.nf — 3′/5′ stranded coverage (+ allMap), per-sample
// ----------------------------------------------------------------------------
// Overview
//   • (Opt) UMI dedup → BAM (+ .bai) used for coverage
//   • Main BAM: 3′ POS/NEG (NEG mirrored); PE also 5′ POS/NEG
//   • allMap BAM: same (3′ always; 5′ if PE)
//   • BigWigs via UCSC bedGraphToBigWig
//   • Parallel per-strand/end; unified step log; README
//
// Why it’s faster
//   • Runs strand/end computations in parallel to use available CPUs
//   • Avoids redundant sorts before BigWig unless forced by a param
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
        umi_tools dedup         -I aligned.bam -S deduplicated.bam --log=${PFX}.dedup_stats.txt
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
  # Clip bedGraph rows to genome bounds and drop unknown contigs
  clip_bg_inplace() {
    local IN_BG="$1" TMP="${1}.cliptmp"
    # Count data rows (non-header) before
    local BEFORE=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "$IN_BG" 2>/dev/null || echo 0)
    awk -v OFS='\t' 'FNR==NR{L[$1]=$2; next}
         (NF>=3)&&($0!~/^(track|browser|#)/){
           chr=$1; s=$2+0; e=$3+0; if(!(chr in L)) next; if(s<0)s=0; m=L[chr]+0; if(m==0) next;
           if(e>m)e=m; if(s<e){ if(NF>=4) print $1,s,e,$4; else print $1,s,e; }
         }' genome.sizes "$IN_BG" > "$TMP" && mv -f "$TMP" "$IN_BG"
    local AFTER=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "$IN_BG" 2>/dev/null || echo 0)
    if [[ "$BEFORE" != "$AFTER" ]]; then
      echo "INFO  [tracks] clip: ${IN_BG} kept=${AFTER} dropped=$(( BEFORE>AFTER ? BEFORE-AFTER : 0 ))"
    else
      echo "INFO  [tracks] clip: ${IN_BG} kept=${AFTER} dropped=0"
    fi
  }

  make_bw () {   # bedGraph → BigWig (clip + optional sort)
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
    
    # POS - run with timeout to prevent hanging
    echo "INFO  [tracks] Computing + strand coverage for ${endflag}' end"
    if ! timeout 1800 bedtools genomecov -bg -strand + -"${endflag}" -ibam "$bam" > "${prefix}.pos.bedgraph"; then
      echo "ERROR  [tracks] Failed to compute + strand coverage (timeout or error)" >&2
      return 1
    fi
    
    # NEG (mirror) - run with timeout
    echo "INFO  [tracks] Computing - strand coverage for ${endflag}' end"
    if ! timeout 1800 bedtools genomecov -bg -strand - -"${endflag}" -ibam "$bam" | awk 'BEGIN{OFS="\\t"}{$4=-$4; print}' > "${prefix}.neg.bedgraph"; then
      echo "ERROR  [tracks] Failed to compute - strand coverage (timeout or error)" >&2
      return 1
    fi
    
    # Clip bedGraph files before BigWig conversion
    echo "INFO  [tracks] Clipping bedGraph files before BigWig conversion"
    clip_bg_inplace "${prefix}.pos.bedgraph"
    clip_bg_inplace "${prefix}.neg.bedgraph"
    
    # Convert to BigWig - run sequentially to avoid wait issues
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

Notes
  - bedtools genomecov with -3 / -5 provides end-specific coverage.
  - Negative-strand tracks are mirrored (score *= -1) for UCSC-style viewing.
  - BigWigs built with bedGraphToBigWig; set params.force_sort_bedgraph=true
    if UCSC sorting complaints appear on your system.
TXT

  echo "INFO  [tracks] ✔ done for ${PFX} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
