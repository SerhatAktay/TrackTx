// ============================================================================
// call_functional_regions.nf — assign 3′ signal to functional classes
// ----------------------------------------------------------------------------
// Overview
//   • Uses siCPM 3′ bedGraphs if present; else CPM.
//   • Active genes only: promoter ∩ divergentTx from input DT BED.
//   • Signal-based counting (bedtools map -o sum) with sequential masking.
//   • Each signal is assigned to exactly ONE functional category hierarchically:
//     Promoter → Gene Body → CPS → Short Genes → Enhancers → Termination Window → Non-localized
//   • Emits BED9 track with per-region counts, a compact summary TSV, README,
//     and a unified step log.
//
// Inputs
//   tuple( sample_id,
//          divergent.bed,
//          pos3_cpm.bg, neg3_cpm.bg,
//          pos3_sicpm.bg, neg3_sicpm.bg,
//          condition, timepoint, replicate )
//   path gtf_file (unused; interface stability)
//   path functional_regions.py
//   path genes.tsv
//   path tss.bed   (optional overrides; BED6 1bp)
//   path tes.bed   (optional overrides; BED6 1bp)
//
// Outputs (publishDir):
//   ${params.output_dir}/07_functional_regions/${sample_id}/
//     ├── functional_regions.bed
//     ├── functional_regions_summary.tsv
//     ├── README_functional_regions.txt
//     └── functional_regions.log
//
// Defaults (params.functional_regions.*):
//   prom_up=250  prom_down=250
//   div_inner=250  div_outer=750
//   tw_length=10000
//   min_signal=0.0
//   allow_unstranded=true
//   count_mode="signal"      // "signal" or "event" (event = legacy 1bp mode)
// ============================================================================

nextflow.enable.dsl = 2

process call_functional_regions {

  // ── Meta / resources ─────────────────────────────────────────────────────
  tag      { sample_id }
  label    'conda'
  cache    'lenient'
  // Resource allocation handled dynamically by base.config

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/07_functional_regions/${sample_id}", mode: 'copy', overwrite: true

  // Use the main tracktx conda environment
  conda (params.conda_fgr ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(divergent_bed),
          path(pos3_cpm_bg), path(neg3_cpm_bg),
          path(pos3_sicpm_bg), path(neg3_sicpm_bg),
          val(condition), val(timepoint), val(replicate)
    path gtf_file
    path functional_regions_py
    path genes_tsv
    path tss_bed
    path tes_bed

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("functional_regions.bed"),
          path("functional_regions_summary.tsv"),
          val(condition), val(timepoint), val(replicate),
          emit: main

    path "README_functional_regions.txt", emit: readme
    path "functional_regions.log",        emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  // Inline all functional_regions defaults directly in the script below to avoid missing-variable issues
  def PROM_UP      = (params.functional_regions?.prom_up      ?: 250)  as int
  def PROM_DOWN    = (params.functional_regions?.prom_down    ?: 250)  as int
  def DIV_INNER    = (params.functional_regions?.div_inner    ?: 250)  as int
  def DIV_OUTER    = (params.functional_regions?.div_outer    ?: 750)  as int
  def TW_LEN       = (params.functional_regions?.tw_length    ?: 10_000) as int
  def MIN_SIGNAL   = (params.functional_regions?.min_signal   ?: 0.0)  as float
  def ALLOW_UNSTR  = (params.functional_regions?.allow_unstranded in [null,true,'true']) ? true : false
  // inline count_mode in script interpolation to avoid scope issues

  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a functional_regions.log) 2>&1
  trap 'echo "ERROR  [FGR] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  

  # ── Bindings / locals ────────────────────────────────────────────────────
  SID="!{sample_id}"
  DIV_BED="!{divergent_bed}"

  POS_CPM="!{pos3_cpm_bg}"; NEG_CPM="!{neg3_cpm_bg}"
  POS_SI="!{pos3_sicpm_bg}"; NEG_SI="!{neg3_sicpm_bg}"

  GENES_TSV="!{genes_tsv}"
  TSS_BED="!{tss_bed}"
  TES_BED="!{tes_bed}"
  FGR_PY="!{functional_regions_py}"

  THREADS=!{task.cpus}

  # Choose 3′ sources (prefer siCPM if non-empty)
  choose_bg(){ local si="$1" cpm="$2"
    [[ -s "$si" ]] && echo "$si" || echo "$cpm"; }
  POS_BG="$(choose_bg "$POS_SI" "$POS_CPM")"
  NEG_BG="$(choose_bg "$NEG_SI" "$NEG_CPM")"

  echo "INFO  [FGR] ▶ sample=\${SID}  cpus=\${THREADS} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "INFO  POS=$(basename "$POS_BG")  NEG=$(basename "$NEG_BG")  count_mode=!{ (params.functional_regions?.count_mode ?: 'signal').toString() }"

  # ── Call Python driver (signal mode default; active genes only) ──────────
  python3 "${FGR_PY}" \
    --sid "${SID}" \
    --genes "${GENES_TSV}" \
    --divergent "${DIV_BED}" \
    --pos "${POS_BG}" \
    --neg "${NEG_BG}" \
    --tss "${TSS_BED}" \
    --tes "${TES_BED}" \
    --prom-up "!{ (params.functional_regions?.prom_up    ?: 250) as int }" \
    --prom-down "!{ (params.functional_regions?.prom_down  ?: 250) as int }" \
    --div-inner "!{ (params.functional_regions?.div_inner  ?: 250) as int }" \
    --div-outer "!{ (params.functional_regions?.div_outer  ?: 750) as int }" \
    --tw-length "!{ (params.functional_regions?.tw_length  ?: 10_000) as int }" \
    --min-signal "!{ (params.functional_regions?.min_signal ?: 0.0) as float }" \
    --count-mode "!{ (params.functional_regions?.count_mode ?: 'signal').toString() }" \
    $([[ "!{ (params.functional_regions?.allow_unstranded in [null,true,'true']) ? true : false }" == "true" ]] && echo "--allow-unstranded" || true) \
    --outdir "." 2>&1 | tee -a functional_regions.log

  # Normalize filenames (driver already writes the canonical names)
  [[ -s functional_regions.bed ]] || : > functional_regions.bed
  [[ -s functional_regions_summary.tsv ]] || printf "region\tsignal\tregion_count\n" > functional_regions_summary.tsv

  # ── README ───────────────────────────────────────────────────────────────
  cat > README_functional_regions.txt <<TXT
Functional-region outputs for sample: ${SID}

Inputs
  - divergent_transcription.bed
  - 3′ tracks: POS=$(basename "$POS_BG"), NEG=$(basename "$NEG_BG")  (count_mode=!{ (params.functional_regions?.count_mode ?: 'signal').toString() })
  - genes.tsv (+ optional TSS/TES overrides)

Geometry (defaults; user-editable via params.functional_regions.*)
  - Promoter:  TSS -!{ (params.functional_regions?.prom_up    ?: 250) as int } .. +!{ (params.functional_regions?.prom_down ?: 250) as int }
  - Divergent: TSS -!{ (params.functional_regions?.div_outer  ?: 750) as int } .. -!{ (params.functional_regions?.div_inner ?: 250) as int }  (opposite strand)
  - CPS:       TES -500 .. +500  (fixed)
  - TW:        CPS end +!{ (params.functional_regions?.tw_length ?: 10_000) as int }
  - Gene body: contiguous between Promoter and CPS

Files
  • functional_regions.bed
  • functional_regions_summary.tsv (region, signal, region_count)
  • functional_regions.log (this step)
TXT

  echo "INFO  [FGR] ✔ done ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
