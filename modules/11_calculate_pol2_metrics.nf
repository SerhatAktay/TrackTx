// ============================================================================
// calculate_pol2_metrics.nf — Pol II metrics (genes+pausing from BAM;
//                              density from normalized 3′ signal)
// ----------------------------------------------------------------------------
// Overview
//   • Density: siCPM (preferred) or CPM 3′ bedGraphs → |pos|+|neg| → sum per
//     functional region via bedtools map. Clear column naming.
//   • Genes & pausing: fast window counts using bedtools coverage on a
//     MAPQ-filtered (and optionally de-duplicated) BAM; no re-sort if not needed.
//   • QC JSON + README; unified step log.
//
// Why it’s faster
//   • Skips resorting when BAM is already SO:coordinate.
//   • Uses bedtools coverage for TSS/body windows (single fast passes).
//
// Inputs
//   tuple( sid, bam, functional.bed|"-",
//          pos3_cpm.bg, neg3_cpm.bg, pos3_sicpm.bg, neg3_sicpm.bg,
//          condition, timepoint, replicate )
//   path gtf
//
// Outputs (publishDir):
//   ${params.output_dir}/08_pol2_metrics/<sid>/
//     ├── pol2_gene_metrics.tsv
//     ├── pausing_index.tsv
//     ├── pol2_density.tsv
//     ├── pol2_qc.json                 (optional)
//     ├── README.txt
//     └── pol2_metrics.log             (unified step log)
// ============================================================================

nextflow.enable.dsl = 2

process calculate_pol2_metrics {

  // ── Meta / resources ─────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  tag    { sid }
  label  'conda'
  cache  'lenient'

  // Use the main tracktx conda environment
  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  publishDir "${params.output_dir}/08_pol2_metrics/${sid}", mode: 'copy', overwrite: true


  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sid),
          path(in_bam),
          path(func_bed),
          path(pos3_cpm_bg), path(neg3_cpm_bg),
          path(pos3_sicpm_bg), path(neg3_sicpm_bg),
          val(cond), val(tp), val(rep)
    path gtf

  // ── Outputs ────────────────────────────────────────────────────────────
  output:
    tuple val(sid), path('pol2_gene_metrics.tsv'), val(cond), val(tp), val(rep), emit: genes
    tuple val(sid), path('pausing_index.tsv'),     val(cond), val(tp), val(rep), emit: pausing
    tuple val(sid), path('pol2_density.tsv'),      val(cond), val(tp), val(rep), emit: density
    path  'pol2_qc.json', optional: true, emit: qc
    path  'README.txt'
    path  'pol2_metrics.log', emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  // Inline params in the script template to avoid scope issues with !{VAR}

  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  trap 'echo "ERROR [pol2] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  exec > >(tee -a pol2_metrics.log) 2>&1
  

  # ── Bindings ─────────────────────────────────────────────────────────────
  SID='!{sid}'
  IN_BAM='!{in_bam}'
  FUNC_BED='!{func_bed}'
  GTF_IN='!{gtf}'
  CALC_PY='!{projectDir}/bin/calculate_pol2_metrics.py'

  POS_CPM='!{pos3_cpm_bg}'; NEG_CPM='!{neg3_cpm_bg}'
  POS_SI='!{pos3_sicpm_bg}'; NEG_SI='!{neg3_sicpm_bg}'

  DEDUP_ON='!{ (params.pol2?.dedup ?: true) ? "true" : "false" }'

  echo "[pol2] ▶ sample=${SID} cpus=!{task.cpus} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  for x in samtools bedtools python3; do
    command -v "$x" >/dev/null 2>&1 || { echo "[pol2] ERROR: missing: $x" >&2; exit 2; }
  done

  [[ -s "$IN_BAM" ]] || { echo "[pol2] ERROR: BAM missing/empty" >&2; exit 2; }
  [[ -e "$GTF_IN" ]] || { echo "[pol2] ERROR: GTF missing" >&2; exit 2; }
  # Require CPM; allow siCPM to be absent (we will fallback per choose_bg)
  for f in "$POS_CPM" "$NEG_CPM"; do [[ -e "$f" ]] || { echo "[pol2] ERROR: bedGraph missing: $f" >&2; exit 2; }; done

  # ── 1) Density from normalized 3′ (siCPM preferred) ──────────────────────
  choose_bg(){
    local si="$1" cpm="$2" out="$3"
    if [[ -s "$si" ]]; then printf "%s" "$si" > "$out"; else printf "%s" "$cpm" > "$out"; fi
  }
  choose_bg "$POS_SI" "$POS_CPM" pos.src
  choose_bg "$NEG_SI" "$NEG_CPM" neg.src
  POS_BG="$(cat pos.src)"; NEG_BG="$(cat neg.src)"
  NORM_METHOD=$( [[ "$POS_BG" == "$POS_SI" ]] && echo siCPM || echo CPM )

  cat_or_zcat(){ case "$1" in *.gz) gzip -cd -- "$1" ;; *) cat -- "$1" ;; esac; }
  abs_clean() {  # in -> out
    local in="$1" out="$2"
    [[ -s "$in" ]] || { : > "$out"; return; }
    cat_or_zcat "$in" | awk 'BEGIN{OFS="\t"} /^[#]|^(track|browser)/{next}
      (NF>=4){ c=$1; s=$2+0; e=$3+0; v=$4+0; if(e>s){ if(v<0)v=-v; print c,s,e,v } }' \
      | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$out"
  }
  abs_clean "$POS_BG" pos.abs.bg
  abs_clean "$NEG_BG" neg.abs.bg

  # Merge to one signal track (|pos|+|neg|)
  cat pos.abs.bg neg.abs.bg \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge -i - -c 4 -o sum > combined.norm.bg || : > combined.norm.bg

  if [[ "$FUNC_BED" != "-" && -s "$FUNC_BED" ]]; then
    CLEAN_FUNC="functional_regions.sorted.bed"
    awk 'BEGIN{OFS="\t"} !/^(track|browser|#)/ && NF>=3 {print $1,$2,$3,(NF>=4?$4:"."),(NF>=5?$5:"0"),(NF>=6?$6:".")}' "$FUNC_BED" \
      | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$CLEAN_FUNC"
    echo -e "chr\tstart\tend\tname\tsignal\tnorm_method" > pol2_density.tsv
    bedtools map -a "$CLEAN_FUNC" -b combined.norm.bg -c 4 -o sum -null 0 \
      | awk -v OFS='\t' -v M="$NORM_METHOD" '{print $1,$2,$3,($4? $4:"."),($5? $5:0),M}' >> pol2_density.tsv
  else
    echo -e "chr\tstart\tend\tname\tsignal\tnorm_method" > pol2_density.tsv
    echo "[pol2] INFO: functional BED absent → emitting density header only."
  fi

  # ── 2) Filter BAM (MAPQ/de-dup) and index; skip re-sort if already SO:coord ─
  SO_COORD=0
  if samtools view -H "$IN_BAM" | awk '/^@HD/ && /SO:coordinate/ {ok=1} END{exit ok?0:1}'; then SO_COORD=1; fi

  DEDUP_FLAG=()
  [[ "$DEDUP_ON" == "true" ]] && DEDUP_FLAG=(-F 0x400)

  if [[ $SO_COORD -eq 1 ]]; then
    echo "[pol2] INFO: BAM appears coordinate-sorted → skip resort"
    samtools view -@ !{task.cpus} -b -q !{ (params.pol2?.mapq ?: 10) as int } -F 0x4 "${DEDUP_FLAG[@]}" "$IN_BAM" -o filtered.bam
  else
    echo "[pol2] INFO: BAM not coordinate-sorted → filtering + sort"
    samtools view -@ !{task.cpus} -b -q !{ (params.pol2?.mapq ?: 10) as int } -F 0x4 "${DEDUP_FLAG[@]}" "$IN_BAM" \
      | samtools sort -@ !{task.cpus} -o filtered.bam
  fi
  samtools index -@ !{task.cpus} filtered.bam

  # ── 3) Pausing & gene metrics (fast coverage path) ───────────────────────
  python3 "${CALC_PY}" \
      --bam filtered.bam \
      --gtf "$GTF_IN" \
      --tss-win !{ (params.pol2?.tss_win ?: 50) as int } \
      --body-offset-min !{ (params.pol2?.body_offset_min ?: 2000) as int } \
      --body-offset-frac !{ (params.pol2?.body_offset_frac ?: 0.10) as float } \
      --feature-types "!{ (params.pol2?.feature_types ?: 'gene,transcript').toString() }" \
      --out-pausing pausing_index.tsv \
      --out-genes   pol2_gene_metrics.tsv \
      --out-qc      pol2_qc.json \
      --threads     "!{task.cpus}" \
      --fail-if-empty "!{ (params.pol2?.fail_if_no_genes ?: false) ? 'true' : 'false' }"

  # Ensure tables exist
  [[ -s pol2_gene_metrics.tsv ]] || echo -e 'gene_id\tgene_name\tchrom\tstrand\ttss_bp\ttss_lo\ttss_hi\ttss_width\tbody_lo\tbody_hi\tbody_len\ttss_count\ttss_cpm\ttss_density_per_bp\tbody_count\tbody_cpm\tbody_density_per_bp\tpi_raw\tpi_len_norm\tis_truncated' > pol2_gene_metrics.tsv
  [[ -s pausing_index.tsv     ]] || echo -e 'gene_id\tchrom\tstrand\tss_count\tgene_body_count\tpausing_index\tis_truncated' > pausing_index.tsv

  # ── 4) README ────────────────────────────────────────────────────────────
  cat > README.txt <<TXT
TrackTx — Pol II metrics (per-sample: ${SID})

Density
  • Normalized 3′ signal: |pos|+|neg| using ${NORM_METHOD}; per-region 'signal' via bedtools map.

Genes & Pausing
  • MAPQ≥!{ (params.pol2?.mapq ?: 10) as int } (${DEDUP_ON:+dedup removed}) from BAM.
  • TSS window ±!{ (params.pol2?.tss_win ?: 50) as int } bp; Body starts at max(!{ (params.pol2?.body_offset_min ?: 2000) as int }, !{ (params.pol2?.body_offset_frac ?: 0.10) as float }×gene_len) to gene end.

Files
  • pol2_gene_metrics.tsv
  • pausing_index.tsv
  • pol2_density.tsv
  • pol2_qc.json (QC summary)
  • pol2_metrics.log (this step)
TXT

  echo "[pol2] ✓ done ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
