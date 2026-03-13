// ============================================================================
// score_enhancer_vs_gene.nf — Enhancer-vs-Gene Scoring for Divergent Sites
// ============================================================================
//
// Purpose:
//   For each sample, take divergent_transcription.bed and gene annotations
//   and compute a continuous enhancer_score ∈ [0,1] per divergent site.
//
// Inputs:
//   tuple(sample_id,
//         divergent_bed,
//         condition, timepoint, replicate)
//   path(genes_tsv)
//   path(tss_bed)
//
// Outputs:
//   ${params.output_dir}/06_divergent_tx/${sample_id}/
//     ├── divergent_enhancer_features.tsv   — Features + rule/ML scores
//     └── divergent_enhancer_score.bed      — BED6: chr, start, end, total, dt_score, enhancer_score
//
// Notes:
//   - Uses bin/enhancer_gene_score.py (Python, numpy, pandas, scikit-learn)
//   - Non-destructive: original divergent_transcription.bed is unchanged.
//   - This process is optional from the standpoint of downstream modules,
//     but the outputs can be surfaced in reports / custom analyses.
//
// ============================================================================

nextflow.enable.dsl = 2

process score_enhancer_vs_gene {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/06_divergent_tx/${sample_id}",
             mode: params.publish_mode,
             overwrite: true

  conda (params.conda_divergent ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(divergent_bed),
          val(condition), val(timepoint), val(replicate)
    path genes_tsv
    path tss_bed

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("divergent_enhancer_score.bed"),
          path("divergent_enhancer_features.tsv"),
          val(condition), val(timepoint), val(replicate),
          emit: scores

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  exec > "enhancer_score.log"
  exec 2> >(tee -a "enhancer_score.log" >&2)

  SAMPLE_ID="!{sample_id}"
  DIV_BED="!{divergent_bed}"
  GENES_TSV="!{genes_tsv}"
  TSS_BED="!{tss_bed}"

  echo "ENHSC | START | sample=${SAMPLE_ID}"
  echo "ENHSC | INPUT | divergent=${DIV_BED}"
  echo "ENHSC | INPUT | genes.tsv=${GENES_TSV}"
  echo "ENHSC | INPUT | tss.bed=${TSS_BED}"

  # Resolve Python interpreter (similar logic as detect_divergent_transcription)
  if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/python3" ]]; then
    PYTHON_CMD="${CONDA_PREFIX}/bin/python3"
  elif [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/python" ]]; then
    PYTHON_CMD="${CONDA_PREFIX}/bin/python"
  elif command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD="python3"
  else
    PYTHON_CMD="python"
  fi

  SCORE_SCRIPT="${projectDir}/bin/enhancer_gene_score.py"
  if [[ ! -f "${SCORE_SCRIPT}" ]]; then
    echo "ENHSC | ERROR | Scoring script not found: ${SCORE_SCRIPT}" >&2
    exit 1
  fi

  # Basic validation: ensure divergent BED exists (may be empty)
  if [[ ! -e "${DIV_BED}" ]]; then
    echo "ENHSC | WARNING | divergent_transcription.bed missing; creating empty outputs"
    : > divergent_enhancer_score.bed
    cat > divergent_enhancer_features.tsv <<EOF
chrom	start	end	center	length	total_signal	dt_score	dist_to_nearest_tss	in_active_tss_window	in_promoter_window	overlaps_gene_body	dist_to_gene_body	upstream_gene_end	downstream_gene_start	intergenic_rel_pos	local_tss_density_100kb	enhancer_score_rule	enhancer_score
EOF
    exit 0
  fi

  OUT_TSV="divergent_enhancer_features.tsv"
  OUT_BED="divergent_enhancer_score.bed"

  echo "ENHSC | RUN   | Computing enhancer-vs-gene scores..."

  set +e
  ${PYTHON_CMD} "${SCORE_SCRIPT}" \
    --divergent "${DIV_BED}" \
    --genes "${GENES_TSV}" \
    --tss "${TSS_BED}" \
    --prom-up "!{params.functional_regions?.prom_up ?: 250}" \
    --prom-down "!{params.functional_regions?.prom_down ?: 250}" \
    --tss-active-pm "!{params.functional_regions?.tss_active_pm ?: 500}" \
    --out-tsv "${OUT_TSV}" \
    --out-bed "${OUT_BED}"
  RC=$?
  set -e

  if [[ ${RC} -ne 0 ]]; then
    echo "ENHSC | ERROR | enhancer_gene_score.py failed with exit code ${RC}" >&2
    exit ${RC}
  fi

  if [[ ! -e "${OUT_BED}" ]]; then
    echo "ENHSC | WARNING | ${OUT_BED} missing after scoring; creating empty file" >&2
    : > "${OUT_BED}"
  fi

  if [[ ! -e "${OUT_TSV}" ]]; then
    echo "ENHSC | WARNING | ${OUT_TSV} missing after scoring; creating stub" >&2
    cat > "${OUT_TSV}" <<EOF
chrom	start	end	center	length	total_signal	dt_score	dist_to_nearest_tss	in_active_tss_window	in_promoter_window	overlaps_gene_body	dist_to_gene_body	upstream_gene_end	downstream_gene_start	intergenic_rel_pos	local_tss_density_100kb	enhancer_score_rule	enhancer_score
EOF
  fi

  echo "ENHSC | COMPLETE | sample=${SAMPLE_ID}"
  '''
}

