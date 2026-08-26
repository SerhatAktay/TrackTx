// ────────────────────────────────────────────────────────────────────────────
// 11b_collect_pol_metrics_per_replicate.nf
// Collect PER-REPLICATE polymerase gene metrics into one tidy long table.
//
// WHY THIS EXISTS
//   When params.replicates.merge=true the main pipeline pools replicates into a
//   single track per condition (n=1). That is correct for visualization, but it
//   destroys the within-condition variance estimate that any differential test
//   needs — so the cohort contrast table collapses to p=1.0 everywhere and is
//   not emitted (see 12_summarize_polymerase_metrics.nf).
//
//   This step instead runs the unchanged gene-metric calculation on each
//   INDIVIDUAL replicate BAM (before merging) and concatenates the results into
//   pol_gene_metrics_per_replicate.tsv. That table — one row per (gene, replicate)
//   — is the intended hand-off for downstream differential analysis in
//   DESeq2 / edgeR, where the user supplies the design and the variance comes
//   from real replication.
//
// INPUT
//   samples_tsv : manifest (sample_id, condition, timepoint, replicate, file)
//                 where `file` is the staged name (metric_N) of that replicate's
//                 pol_gene_metrics.tsv
//   metric_*    : the staged per-replicate pol_gene_metrics.tsv files
//
// OUTPUT
//   pol_gene_metrics_per_replicate.tsv  (long format; columns below)
// ────────────────────────────────────────────────────────────────────────────

process collect_pol_metrics_per_replicate {

  tag        'per-replicate'
  label      'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/08b_pol_metrics_per_replicate" },
             mode: params.publish_mode,
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  input:
    path(samples_tsv, stageAs: 'samples.tsv')
    path('metric_*')

  output:
    path 'pol_gene_metrics_per_replicate.tsv', emit: table
    path 'README_per_replicate.txt'
    path 'collect_per_replicate.log',          emit: log

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  exec > >(tee -a collect_per_replicate.log)
  exec 2> >(tee -a collect_per_replicate.log >&2)

  echo "PER-REPLICATE | START | Collecting per-replicate Pol-II gene metrics"

  OUT='pol_gene_metrics_per_replicate.tsv'
  printf 'gene_id\\tgene_name\\tsample_id\\tcondition\\ttimepoint\\treplicate\\ttss_cpm\\tbody_cpm\\tpi_raw\\tpi_len_norm\\n' > "\${OUT}"

  # Nextflow stages a single file as 'metric_' (no index). Normalize so the
  # manifest's metric_N references always resolve.
  if [[ -e metric_ && ! -e metric_1 ]]; then mv metric_ metric_1; fi

  # Skip the manifest header; iterate one replicate per line.
  tail -n +2 samples.tsv | while IFS=\$'\\t' read -r SID COND TP REP FNAME; do
    [[ -z "\${SID:-}" ]] && continue
    if [[ ! -s "\${FNAME}" ]]; then
      echo "PER-REPLICATE | WARN | Missing/empty metric file for \${SID} (\${FNAME}) — skipping"
      continue
    fi
    # Pull the four metric columns BY HEADER NAME (column order is not assumed),
    # prepend sample metadata, append to the long table.
    awk -F'\\t' -v sid="\${SID}" -v cond="\${COND}" -v tp="\${TP}" -v rep="\${REP}" '
      NR==1 {
        for (i=1;i<=NF;i++) h[\$i]=i
        next
      }
      {
        gi  = (("gene_id"      in h) ? \$h["gene_id"]      : "")
        gn  = (("gene_name"    in h) ? \$h["gene_name"]    : "")
        tc  = (("tss_cpm"      in h) ? \$h["tss_cpm"]      : "")
        bc  = (("body_cpm"     in h) ? \$h["body_cpm"]     : "")
        pr  = (("pi_raw"       in h) ? \$h["pi_raw"]       : "")
        pln = (("pi_len_norm"  in h) ? \$h["pi_len_norm"]  : "")
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", gi, gn, sid, cond, tp, rep, tc, bc, pr, pln
      }
    ' "\${FNAME}" >> "\${OUT}"
  done

  N_LINES=\$(( \$(wc -l < "\${OUT}") - 1 ))
  echo "PER-REPLICATE | DONE | Wrote \${N_LINES} gene-rows to \${OUT}"

  cat > README_per_replicate.txt <<'DOCEOF'
================================================================================
PER-REPLICATE POLYMERASE GENE METRICS — pol_gene_metrics_per_replicate.tsv
================================================================================

PURPOSE
  One row per (gene, replicate). This is the analytical hand-off for differential
  testing. The merged cohort outputs (09_pol_aggregate/) pool replicates to n=1
  per condition and therefore cannot support per-gene statistics; use THIS table
  instead when you want p-values.

COLUMNS
  gene_id       Ensembl/RefSeq gene id
  gene_name     gene symbol
  sample_id     individual replicate sample id (pre-merge)
  condition     experimental condition
  timepoint     timepoint label
  replicate     replicate number
  tss_cpm       TSS-window signal, CPM
  body_cpm      gene-body signal, CPM
  pi_raw        raw pausing index (TSS density / body density)
  pi_len_norm   length-normalized pausing index

SUGGESTED DOWNSTREAM USE (R)
  library(tidyverse)
  d <- read_tsv("pol_gene_metrics_per_replicate.tsv")
  # e.g. build a counts-like matrix from body_cpm and run a paired comparison in
  # DESeq2/edgeR with ~ condition (+ replicate), with >=2 replicates per group.

NOTE
  Metrics are computed with the SAME calculate_pol_metrics.py call used for the
  merged tracks, run on each individual replicate BAM, so values are directly
  comparable to 08_pol_metrics/ (which are the merged equivalents).
================================================================================
DOCEOF

  echo "PER-REPLICATE | END"
  """
}
