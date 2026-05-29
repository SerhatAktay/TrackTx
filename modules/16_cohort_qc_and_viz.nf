// ============================================================================
// cohort_qc_and_viz.nf — MultiQC, deepTools QC, IGV Session, Run-on Efficiency
// ============================================================================
//
// Purpose:
//   Final cohort-level QC and visualization module that runs once all samples
//   are complete.  Produces:
//     1. MultiQC HTML report aggregating all per-sample QC logs
//     2. deepTools PCA plot and Pearson-correlation heatmap (BigWig-based)
//     3. IGV session XML for one-click track loading
//     4. Run-on efficiency table (5'/3' signal ratio per gene body)
//
// Run-on efficiency:
//   PRO-seq run-on produces signal across the full gene body.  Poor run-on
//   shows a 5'/3' signal drop — reads pile up near the TSS and fall off.
//   Metric: median(5p_body_signal / 3p_body_signal) per sample, computed
//   over a set of long, expressed genes (>= 10 kb, 3p signal > threshold).
//   Values close to 1.0 indicate efficient run-on; <0.3 suggests poor run-on.
//
// deepTools PCA / correlation:
//   Operates on CPM-normalized 3' positive-strand BigWigs (one per sample).
//   Uses 10 kb genome bins (default) for fast genome-wide summary.
//
// IGV session:
//   Groups tracks by sample, organises as:
//     <sample>  3p pos (CPM) / 3p neg (CPM) / allMap pos / allMap neg
//   Colour-codes conditions automatically (up to 8 distinct colours).
//
// Inputs:
//   path(multiqc_log_dir) : directory tree containing all QC logs
//   val(bw_pos3_list)     : list of CPM 3' positive BigWig paths (all samples)
//   val(bw_neg3_list)     : list of CPM 3' negative BigWig paths (all samples)
//   val(bw_allmap_pos_list): allMap positive BigWig paths
//   val(bw_allmap_neg_list): allMap negative BigWig paths
//   val(sample_ids)        : sample identifiers (same order as BW lists)
//   val(conditions)        : condition labels (same order)
//   val(pos3_bg_list)      : raw 3' positive bedGraphs for run-on calculation
//   val(neg3_bg_list)      : raw 3' negative bedGraphs
//   val(pos5_bg_list)      : raw 5' positive bedGraphs
//   val(neg5_bg_list)      : raw 5' negative bedGraphs
//   path(genes_bed)        : BED6 gene annotation (from download_genome_annotations)
//
// Outputs:
//   ${params.output_dir}/12_cohort_qc/
//     ├── multiqc/
//     │   ├── multiqc_report.html
//     │   └── multiqc_data/
//     ├── deeptools/
//     │   ├── bigwig_summary.npz
//     │   ├── pca_plot.pdf
//     │   └── correlation_heatmap.pdf
//     ├── igv_session.xml
//     └── runon_efficiency.tsv
//
// ============================================================================


process cohort_qc_and_viz {

  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/12_cohort_qc",
             mode: params.publish_mode,
             overwrite: true

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    path  multiqc_log_dir         // flat collection of all QC log files
    val   bw_pos3_list
    val   bw_neg3_list
    val   bw_allmap_pos_list
    val   bw_allmap_neg_list
    val   sample_ids
    val   conditions
    val   pos3_bg_list
    val   neg3_bg_list
    val   pos5_bg_list
    val   neg5_bg_list
    path  genes_bed

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    path "multiqc/multiqc_report.html",            emit: multiqc_html,    optional: true
    path "multiqc/multiqc_data",                   emit: multiqc_data,    optional: true
    path "deeptools/bigwig_summary.npz",           emit: bw_summary,      optional: true
    path "deeptools/pca_plot.pdf",                 emit: pca_plot,        optional: true
    path "deeptools/correlation_heatmap.pdf",      emit: corr_heatmap,    optional: true
    path "igv_session.xml",                        emit: igv_session
    path "runon_efficiency.tsv",                   emit: runon_efficiency

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  // Groovy: materialise val lists into shell-accessible strings
  bwPos3Str      = (bw_pos3_list      instanceof List ? bw_pos3_list      : [bw_pos3_list])
                     .findAll { it && it != 'null' }.join(' ')
  bwNeg3Str      = (bw_neg3_list      instanceof List ? bw_neg3_list      : [bw_neg3_list])
                     .findAll { it && it != 'null' }.join(' ')
  bwAmPos3Str    = (bw_allmap_pos_list instanceof List ? bw_allmap_pos_list : [bw_allmap_pos_list])
                     .findAll { it && it != 'null' }.join(' ')
  bwAmNeg3Str    = (bw_allmap_neg_list instanceof List ? bw_allmap_neg_list : [bw_allmap_neg_list])
                     .findAll { it && it != 'null' }.join(' ')
  sampleIdsStr   = (sample_ids  instanceof List ? sample_ids  : [sample_ids])
                     .findAll { it && it != 'null' }.join(' ')
  conditionsStr  = (conditions  instanceof List ? conditions  : [conditions])
                     .findAll { it && it != 'null' }.join(' ')
  pos3BgStr      = (pos3_bg_list instanceof List ? pos3_bg_list : [pos3_bg_list])
                     .findAll { it && it != 'null' }.join(' ')
  neg3BgStr      = (neg3_bg_list instanceof List ? neg3_bg_list : [neg3_bg_list])
                     .findAll { it && it != 'null' }.join(' ')
  pos5BgStr      = (pos5_bg_list instanceof List ? pos5_bg_list : [pos5_bg_list])
                     .findAll { it && it != 'null' }.join(' ')
  neg5BgStr      = (neg5_bg_list instanceof List ? neg5_bg_list : [neg5_bg_list])
                     .findAll { it && it != 'null' }.join(' ')
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  exec > >(tee -a cohort_qc.log)
  exec 2> >(tee -a cohort_qc.log >&2)

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COHORT_QC | START | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  # Materialise shell arrays from Nextflow-interpolated strings
  read -ra BW_POS3       <<< "!{bwPos3Str}"
  read -ra BW_NEG3       <<< "!{bwNeg3Str}"
  read -ra BW_AMPOS3     <<< "!{bwAmPos3Str}"
  read -ra BW_AMNEG3     <<< "!{bwAmNeg3Str}"
  read -ra SAMPLE_IDS    <<< "!{sampleIdsStr}"
  read -ra CONDITIONS    <<< "!{conditionsStr}"
  read -ra POS3_BGS      <<< "!{pos3BgStr}"
  read -ra NEG3_BGS      <<< "!{neg3BgStr}"
  read -ra POS5_BGS      <<< "!{pos5BgStr}"
  read -ra NEG5_BGS      <<< "!{neg5BgStr}"
  GENES_BED="!{genes_bed}"
  N_SAMPLES="${#SAMPLE_IDS[@]}"

  echo "COHORT_QC | CONFIG | Samples: ${N_SAMPLES}"
  echo "COHORT_QC | CONFIG | Sample IDs: ${SAMPLE_IDS[*]}"

  mkdir -p multiqc deeptools

  ###########################################################################
  # 1) MULTIQC
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COHORT_QC | MULTIQC | Aggregating QC logs..."
  echo "────────────────────────────────────────────────────────────────────────"

  if command -v multiqc >/dev/null 2>&1; then
    # multiqc_log_dir is a flat collection of staged log files
    multiqc . \
      --outdir multiqc \
      --filename multiqc_report.html \
      --title "TrackTx Cohort QC" \
      --quiet \
      2>&1 | sed 's/^/COHORT_QC | MULTIQC | /'
    echo "COHORT_QC | MULTIQC | Report created: multiqc/multiqc_report.html"
  else
    echo "COHORT_QC | MULTIQC | WARNING: multiqc not found — skipping"
    echo "COHORT_QC | MULTIQC | Install via: pip install multiqc"
  fi

  ###########################################################################
  # 2) DEEPTOOLS: PCA + CORRELATION HEATMAP
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COHORT_QC | DEEPTOOLS | Building BigWig summary for PCA/correlation..."
  echo "────────────────────────────────────────────────────────────────────────"

  DT_OK=0
  if command -v multiBigwigSummary >/dev/null 2>&1 && [[ ${N_SAMPLES} -ge 2 ]]; then
    # Filter to BigWigs that actually exist on disk
    VALID_BW=()
    for bw in "${BW_POS3[@]}"; do
      [[ -s "${bw}" ]] && VALID_BW+=("${bw}")
    done

    if [[ ${#VALID_BW[@]} -ge 2 ]]; then
      echo "COHORT_QC | DEEPTOOLS | Using ${#VALID_BW[@]} BigWig files"

      multiBigwigSummary bins \
        -b "${VALID_BW[@]}" \
        -l "${SAMPLE_IDS[@]}" \
        -o deeptools/bigwig_summary.npz \
        --binSize 10000 \
        --numberOfProcessors !{task.cpus} \
        --outRawCounts deeptools/bigwig_summary_counts.tsv \
        --skipZeroOverZero \
        2>&1 | sed 's/^/COHORT_QC | DEEPTOOLS | /' && DT_OK=1

      if [[ ${DT_OK} -eq 1 ]]; then
        echo "COHORT_QC | DEEPTOOLS | Generating PCA plot..."
        plotPCA \
          -in deeptools/bigwig_summary.npz \
          -o  deeptools/pca_plot.pdf \
          --plotTitle "TrackTx PCA — 3' CPM signal" \
          --outFileNameData deeptools/pca_data.tsv \
          2>&1 | sed 's/^/COHORT_QC | DEEPTOOLS | PCA | /' || true

        echo "COHORT_QC | DEEPTOOLS | Generating correlation heatmap..."
        plotCorrelation \
          -in deeptools/bigwig_summary.npz \
          --corMethod pearson \
          --skipZeros \
          --whatToPlot heatmap \
          --colorMap RdYlBu_r \
          --plotNumbers \
          -o deeptools/correlation_heatmap.pdf \
          --outFileCorMatrix deeptools/correlation_matrix.tsv \
          2>&1 | sed 's/^/COHORT_QC | DEEPTOOLS | CORR | /' || true

        echo "COHORT_QC | DEEPTOOLS | Complete"
      fi
    else
      echo "COHORT_QC | DEEPTOOLS | Fewer than 2 valid BigWigs found — skipping"
    fi
  else
    if ! command -v multiBigwigSummary >/dev/null 2>&1; then
      echo "COHORT_QC | DEEPTOOLS | WARNING: deeptools not found — skipping"
    else
      echo "COHORT_QC | DEEPTOOLS | Only 1 sample — PCA/correlation requires >=2"
    fi
  fi

  # Ensure placeholder outputs exist so publishDir does not fail
  [[ -s deeptools/bigwig_summary.npz ]]       || touch deeptools/bigwig_summary.npz
  [[ -s deeptools/pca_plot.pdf ]]             || touch deeptools/pca_plot.pdf
  [[ -s deeptools/correlation_heatmap.pdf ]]  || touch deeptools/correlation_heatmap.pdf

  ###########################################################################
  # 3) IGV SESSION XML
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COHORT_QC | IGV | Generating IGV session file..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Colour palette: up to 8 distinct conditions get distinct colours.
  # Tracks for the same condition share a colour.
  PALETTE=("31,119,180" "255,127,14" "44,160,44" "214,39,40"
           "148,103,189" "140,86,75" "227,119,194" "127,127,127")

  python3 - \
    "${SAMPLE_IDS[@]}" \
    "---conditions" \
    "${CONDITIONS[@]}" \
    "---bw_pos3" \
    "${BW_POS3[@]}" \
    "---bw_neg3" \
    "${BW_NEG3[@]}" \
    "---bw_ampos3" \
    "${BW_AMPOS3[@]}" \
    "---bw_amneg3" \
    "${BW_AMNEG3[@]}" \
    <<'PYEOF'
import sys, xml.dom.minidom

args = sys.argv[1:]

# Parse sentinel-separated sections
def split_by(lst, sentinel):
    result, cur = [], []
    for x in lst:
        if x == sentinel:
            result.append(cur); cur = []
        else:
            cur.append(x)
    result.append(cur)
    return result

sections = split_by(args, '---conditions')
sample_ids   = sections[0]
rest         = sections[1] if len(sections) > 1 else []
rest2        = split_by(rest, '---bw_pos3')
conditions   = rest2[0]
rest3        = split_by(rest2[1] if len(rest2) > 1 else [], '---bw_neg3')
bw_pos3      = rest3[0]
rest4        = split_by(rest3[1] if len(rest3) > 1 else [], '---bw_ampos3')
bw_neg3      = rest4[0]
rest5        = split_by(rest4[1] if len(rest4) > 1 else [], '---bw_amneg3')
bw_ampos3    = rest5[0]
bw_amneg3    = rest5[1] if len(rest5) > 1 else []

palette = [
    "31,119,180", "255,127,14", "44,160,44", "214,39,40",
    "148,103,189", "140,86,75", "227,119,194", "127,127,127"
]

unique_conds = []
for c in conditions:
    if c not in unique_conds:
        unique_conds.append(c)
cond_color = {c: palette[i % len(palette)] for i, c in enumerate(unique_conds)}

lines = ['<?xml version="1.0" encoding="UTF-8"?>',
         '<Session genome="hg38" version="8">',
         '  <Resources>']

def safe_path(p):
    """Return absolute or relative path usable in IGV."""
    if p and p != 'null':
        return p
    return None

for i, sid in enumerate(sample_ids):
    cond  = conditions[i] if i < len(conditions) else "unknown"
    color = cond_color.get(cond, "127,127,127")

    p3    = safe_path(bw_pos3[i]   if i < len(bw_pos3)   else None)
    n3    = safe_path(bw_neg3[i]   if i < len(bw_neg3)   else None)
    ap3   = safe_path(bw_ampos3[i] if i < len(bw_ampos3) else None)
    an3   = safe_path(bw_amneg3[i] if i < len(bw_amneg3) else None)

    for path, label in [(p3,  f"{sid} | 3p pos (CPM)"),
                        (n3,  f"{sid} | 3p neg (CPM)"),
                        (ap3, f"{sid} | allMap 3p pos"),
                        (an3, f"{sid} | allMap 3p neg")]:
        if path:
            lines.append(f'    <Resource path="{path}" label="{label}" color="{color}"/>')

lines += ['  </Resources>', '  <Panel name="DataPanel">']

for i, sid in enumerate(sample_ids):
    cond  = conditions[i] if i < len(conditions) else "unknown"
    color = cond_color.get(cond, "127,127,127")
    p3    = safe_path(bw_pos3[i]   if i < len(bw_pos3)   else None)
    n3    = safe_path(bw_neg3[i]   if i < len(bw_neg3)   else None)
    ap3   = safe_path(bw_ampos3[i] if i < len(bw_ampos3) else None)
    an3   = safe_path(bw_amneg3[i] if i < len(bw_amneg3) else None)

    for path, label, yscale in [
        (p3,  f"{sid} | 3p pos (CPM)",     "0,10"),
        (n3,  f"{sid} | 3p neg (CPM)",    "-10,0"),
        (ap3, f"{sid} | allMap 3p pos",    "0,10"),
        (an3, f"{sid} | allMap 3p neg",   "-10,0"),
    ]:
        if path:
            lo, hi = yscale.split(',')
            lines.append(
                f'    <Track id="{path}" name="{label}" color="{color}" '
                f'renderer="BAR_CHART" height="50" '
                f'dataRange="{lo},{hi}" autoscale="false"/>'
            )

lines += ['  </Panel>', '</Session>']

with open('igv_session.xml', 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"IGV session written for {len(sample_ids)} sample(s), {len(unique_conds)} condition(s)")
PYEOF

  echo "COHORT_QC | IGV | igv_session.xml created"

  ###########################################################################
  # 4) RUN-ON EFFICIENCY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COHORT_QC | RUNON | Calculating run-on efficiency..."
  echo "────────────────────────────────────────────────────────────────────────"

  # Metric: for each sample, compute median(5p_body / 3p_body) over long genes.
  # Uses raw (not CPM-normalised) bedGraphs so the ratio cancels out library
  # size — the ratio itself is the diagnostic.
  #
  # Gene body window: TSS+500 → TES-500 (avoids promoter pile-up at 5').
  # Genes included: length >= 10,000 bp, BED6 col 6 = strand.
  #
  # For the negative strand: bedgraph values are negative (mirrored), so we
  # take abs() before summing.

  python3 - \
    "${GENES_BED}" \
    "${#SAMPLE_IDS[@]}" \
    "${SAMPLE_IDS[@]}" \
    "${POS3_BGS[@]}" \
    "${NEG3_BGS[@]}" \
    "${POS5_BGS[@]}" \
    "${NEG5_BGS[@]}" \
    <<'PYEOF'
import sys
import os

args      = sys.argv[1:]
genes_bed = args[0]
n         = int(args[1])
sids      = args[2:2+n]
pos3      = args[2+n  :2+2*n]
neg3      = args[2+2*n:2+3*n]
pos5      = args[2+3*n:2+4*n]
neg5      = args[2+4*n:2+5*n]

MIN_GENE_LEN = 10_000
TSS_SKIP     = 500    # skip first 500 bp after TSS
TES_SKIP     = 500    # skip last  500 bp before TES
MIN_SIGNAL   = 1.0    # minimum 3' signal to include gene

def load_bedgraph(path):
    """Return dict: chrom -> sorted list of (start, end, abs_val)."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return {}
    chrom_data = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom, start, end, val = parts[0], int(parts[1]), int(parts[2]), abs(float(parts[3]))
            chrom_data.setdefault(chrom, []).append((start, end, val))
    for lst in chrom_data.values():
        lst.sort()
    return chrom_data

def region_signal(bg, chrom, lo, hi):
    """Sum bedGraph signal in [lo, hi)."""
    if chrom not in bg:
        return 0.0
    total = 0.0
    for (s, e, v) in bg[chrom]:
        if e <= lo:
            continue
        if s >= hi:
            break
        overlap = min(e, hi) - max(s, lo)
        if overlap > 0:
            total += v * overlap
    return total

# Parse gene bodies from BED6
genes = []
if os.path.isfile(genes_bed):
    with open(genes_bed) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 6:
                continue
            chrom, start, end, name, _, strand = p[0], int(p[1]), int(p[2]), p[3], p[4], p[5]
            length = end - start
            if length < MIN_GENE_LEN:
                continue
            if strand == '+':
                body_lo = start + TSS_SKIP
                body_hi = end   - TES_SKIP
            else:
                body_lo = start + TES_SKIP
                body_hi = end   - TSS_SKIP
            if body_hi <= body_lo:
                continue
            genes.append((chrom, body_lo, body_hi, strand))
else:
    print(f"WARNING: genes_bed not found: {genes_bed}", file=sys.stderr)

print(f"Loaded {len(genes)} gene bodies (>={MIN_GENE_LEN} bp)", file=sys.stderr)

header = ['sample_id', 'n_genes_used', 'median_5p3p_ratio',
          'mean_5p3p_ratio', 'interpretation']
rows = ['\t'.join(header)]

for i, sid in enumerate(sids):
    p3bg = load_bedgraph(pos3[i] if i < len(pos3) else '')
    n3bg = load_bedgraph(neg3[i] if i < len(neg3) else '')
    p5bg = load_bedgraph(pos5[i] if i < len(pos5) else '')
    n5bg = load_bedgraph(neg5[i] if i < len(neg5) else '')

    ratios = []
    for (chrom, lo, hi, strand) in genes:
        if strand == '+':
            sig3 = region_signal(p3bg, chrom, lo, hi)
            sig5 = region_signal(p5bg, chrom, lo, hi)
        else:
            sig3 = region_signal(n3bg, chrom, lo, hi)
            sig5 = region_signal(n5bg, chrom, lo, hi)
        if sig3 < MIN_SIGNAL:
            continue
        ratios.append(sig5 / sig3 if sig3 > 0 else 0.0)

    if ratios:
        ratios.sort()
        mid = len(ratios) // 2
        med = ratios[mid] if len(ratios) % 2 == 1 else (ratios[mid-1] + ratios[mid]) / 2
        avg = sum(ratios) / len(ratios)
        if   med >= 0.7: interp = "excellent"
        elif med >= 0.4: interp = "good"
        elif med >= 0.2: interp = "moderate — check run-on time"
        else:            interp = "poor — consider longer run-on"
    else:
        med, avg, interp = float('nan'), float('nan'), "insufficient_data"

    rows.append(f"{sid}\t{len(ratios)}\t{med:.4f}\t{avg:.4f}\t{interp}")
    print(f"  {sid}: median 5p/3p = {med:.4f} ({interp})", file=sys.stderr)

with open('runon_efficiency.tsv', 'w') as f:
    f.write('\n'.join(rows) + '\n')

print("Run-on efficiency table written.", file=sys.stderr)
PYEOF

  echo "COHORT_QC | RUNON | runon_efficiency.tsv created"

  ###########################################################################
  # FINAL SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COHORT_QC | SUMMARY"
  echo "────────────────────────────────────────────────────────────────────────"
  [[ -s multiqc/multiqc_report.html ]] \
    && echo "COHORT_QC | SUMMARY | ✓ MultiQC report" \
    || echo "COHORT_QC | SUMMARY | ✗ MultiQC report (skipped)"
  [[ -s deeptools/pca_plot.pdf ]] \
    && echo "COHORT_QC | SUMMARY | ✓ deepTools PCA" \
    || echo "COHORT_QC | SUMMARY | ✗ deepTools PCA (skipped)"
  [[ -s deeptools/correlation_heatmap.pdf ]] \
    && echo "COHORT_QC | SUMMARY | ✓ deepTools correlation heatmap" \
    || echo "COHORT_QC | SUMMARY | ✗ deepTools correlation heatmap (skipped)"
  echo "COHORT_QC | SUMMARY | ✓ IGV session XML"
  echo "COHORT_QC | SUMMARY | ✓ Run-on efficiency table"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COHORT_QC | COMPLETE | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}
