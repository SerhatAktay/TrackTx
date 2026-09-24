// ============================================================================
// normalize_coverage_tracks.nf — Track Normalization (CPM and Spike-in CPM)
// ============================================================================
//
// Purpose:
//   Normalizes raw coverage tracks to CPM and spike-in CPM (siCPM)
//
// Features:
//   • CPM normalization: Counts per million mapped reads
//   • siCPM normalization: Spike-in normalized CPM
//   • Processes both main and allMap BAM tracks
//   • Handles 3' tracks (always) and 5' tracks (always)
//   • Single-pass scaling for efficiency
//   • Optional BigWig generation
//   • Comprehensive manifest of all outputs
//
// Normalization Methods:
//   CPM   = (count / sample_reads) × 1,000,000
//   siCPM = (count / sample_spike) × (control_spike / control_reads) × 1,000,000
//
// siCPM Control Selection:
//   1) Condition matches params.control_label (case-insensitive), lowest
//      replicate number preferred (so a merged rep-0 track wins over stray
//      per-replicate rows). Set params.control_label explicitly per dataset
//      to override auto-detection or to disambiguate a multi-arm design
//      (see step 2).
//   2) AUTO-DETECT (used when control_label doesn't match anything): every
//      TrackTx samplesheet has a structured timepoint column, and the
//      baseline/untreated sample is the minimum timepoint in the
//      timecourse (0 in every dataset shipped with this pipeline) —
//      independent of whatever free-text condition name a given study
//      uses. This is trusted ONLY when unambiguous: exactly one distinct
//      condition value at that minimum timepoint with spike reads. A
//      design with more than one arm starting at timepoint 0 (e.g. a
//      "primed" vs "unprimed" pre-treatment) is genuinely ambiguous —
//      auto-detection refuses rather than guessing which arm a given
//      sample belongs to, and control_label must be set by hand for that
//      dataset.
//   If neither resolves a control, siCPM is DISABLED (factor = 0) with a
//   loud, actionable warning — there is deliberately no "guess a control
//   from row order" fallback. An earlier version of this pipeline silently
//   fell back to "first sample with spike_reads > 0" when the label didn't
//   match, which produced a real but wrong siCPM factor (using an
//   arbitrary, often treated-not-baseline sample as the reference) without
//   any indication anything was wrong — see git history / project notes
//   for the incident that surfaced this. Failing safe to "disabled and
//   visibly zero" is far preferable to a plausible-looking but incorrect
//   number; auto-detection only fires when it can be exactly as certain as
//   an explicit label match.
//
// Inputs:
//   tuple(sample_id, pos3_bg, neg3_bg, pos5_bg, neg5_bg,
//         allmap3p_pos_bg, allmap3p_neg_bg, allmap5p_pos_bg, allmap5p_neg_bg,
//         condition, timepoint, replicate, counts_master_tsv, genes_tsv)
//   path(genome_fa) : Genome FASTA (for chromosome sizes)
//
// Outputs:
//   ${params.output_dir}/05_normalized_tracks/${sample_id}/
//     ├── 3p/
//     │   ├── ${sample_id}.3p.pos.cpm.bedgraph
//     │   ├── ${sample_id}.3p.neg.cpm.bedgraph
//     │   ├── ${sample_id}.3p.pos.sicpm.bedgraph
//     │   ├── ${sample_id}.3p.neg.sicpm.bedgraph
//     │   ├── ${sample_id}.3p.pos.cpm.bw
//     │   ├── ${sample_id}.3p.neg.cpm.bw
//     │   ├── ${sample_id}.3p.pos.sicpm.bw
//     │   ├── ${sample_id}.3p.neg.sicpm.bw
//     │   ├── ${sample_id}.allMap.3p.*.bedgraph (if emit_allmap)
//     │   └── ${sample_id}.allMap.3p.*.bw (if emit_allmap)
//     ├── 5p/ (if emit_5p or auto-detected)
//     │   └── [same structure as 3p/]
//     ├── normalization_factors.tsv
//     ├── tracks_manifest.tsv
//     ├── README_normalization.txt
//     └── normalize_coverage_tracks.log
//
// Parameters:
//   params.norm.emit_bw       : Generate BigWig files (default: true)
//   params.norm.emit_sicpm    : Generate siCPM tracks (default: true)
//   params.norm.emit_allmap   : Include allMap tracks (default: true)
//   params.norm.emit_5p       : Force 5' track generation (default: auto)
//   params.force_sort_bedgraph: Sort before BigWig (default: false)
//   params.control_label      : Control condition name — if set, MUST match
//                                a real condition in your samplesheet (e.g.
//                                "no_heat_shock"); default "CTRL" is a
//                                placeholder. When it doesn't match, the
//                                minimum-timepoint sample is used instead
//                                IF it is unambiguous (see siCPM Control
//                                Selection above) — set control_label
//                                explicitly for any multi-arm design.
//   params.norm.timeout_bw    : BigWig timeout seconds (default: 900)
//   params.norm.gene_end_method   : none | tes_window | gene_body (default: none)
//                                    Alternative CPM denominator computed from
//                                    raw 3' signal in a per-gene end region
//                                    (written to normalization_factors.tsv as
//                                    gene_end_<method>; no track is scaled by it).
//   params.norm.gene_end_window   : bp window for tes_window method (default: 500)
//   params.norm.gene_end_min_reads: min reads/gene to count toward the total (default: 10)
//
// ============================================================================


process normalize_coverage_tracks {

  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/05_normalized_tracks/${sample_id}" },
             mode: params.publish_mode,
             overwrite: true,
             saveAs: { filename ->
               def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
               // Skip bedGraphs when output.bedgraph: false (BigWigs sufficient for genome browsers)
               if (params.get('output')?.get('bedgraph')?.toString() == 'false' && name.endsWith('.bedgraph')) return null
               return name
             }

  conda (params.conda_norm ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  // NOTE: bedGraph inputs MUST be path() not val() for proper Docker staging
  input:
    tuple val(sample_id),
          path(pos3_bg),
          path(neg3_bg),
          path(pos5_bg),
          path(neg5_bg),
          path(am3p_pos_bg),
          path(am3p_neg_bg),
          path(am5p_pos_bg),
          path(am5p_neg_bg),
          val(condition), val(timepoint), val(replicate),
          path(counts_master_tsv),
          path(genes_tsv)
    path genome_fa

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    // Legacy compatibility tuple
    tuple val(sample_id),
          path("3p/${sample_id}_pos3_cpm.bedgraph"),
          path("3p/${sample_id}_neg3_cpm.bedgraph"),
          path("normalization_factors.tsv"),
          val(condition), val(timepoint), val(replicate),
          emit: norm_tuple

    // Main CPM tracks
    path "3p/${sample_id}.3p.pos.cpm.bedgraph", emit: pos3_cpm_bg
    path "3p/${sample_id}.3p.neg.cpm.bedgraph", emit: neg3_cpm_bg

    // siCPM 3' bedGraphs, keyed by sample for joining downstream. These files are
    // ALWAYS produced (real siCPM values, or an empty file when siCPM is
    // unavailable), so they are safe to stage. Emitting them here lets module 11
    // receive siCPM via the Nextflow channel (work dir) instead of reading the
    // publish dir — which was both pointing at a wrong path AND would break when
    // output.bedgraph=false (bedGraphs not published).
    tuple val(sample_id),
          path("3p/${sample_id}.3p.pos.sicpm.bedgraph"),
          path("3p/${sample_id}.3p.neg.sicpm.bedgraph"),
          emit: sicpm3p_bg

    // Main 5' CPM/siCPM bedGraphs, same "always produced (real values, or an
    // empty file)" contract as the 3' emits above. GRO-seq's Pol II position
    // is the 5' end (see nextflow.config's signal_end doc), so STEP 11/12 in
    // main.nf need these on a channel to route functional-region/Pol-density
    // signal through the correct end instead of always using 3'.
    tuple val(sample_id),
          path("5p/${sample_id}.5p.pos.cpm.bedgraph"),
          path("5p/${sample_id}.5p.neg.cpm.bedgraph"),
          emit: cpm5p_bg

    tuple val(sample_id),
          path("5p/${sample_id}.5p.pos.sicpm.bedgraph"),
          path("5p/${sample_id}.5p.neg.sicpm.bedgraph"),
          emit: sicpm5p_bg

    // CPM 3' BigWigs handed to the per-sample report (module 14) via channel, so
    // its track-link availability no longer depends on probing the publish dir.
    // All four ALWAYS exist (0-byte placeholder when a track wasn't produced;
    // the report treats 0-byte as "not available"). See section 9c.
    tuple val(sample_id),
          path("3p/${sample_id}.3p.pos.cpm.bw"),
          path("3p/${sample_id}.3p.neg.cpm.bw"),
          path("3p/${sample_id}.allMap.3p.pos.cpm.bw"),
          path("3p/${sample_id}.allMap.3p.neg.cpm.bw"),
          emit: report_bw

    // Documentation and manifest
    path "README_normalization.txt", emit: readme
    path "tracks_manifest.tsv",      emit: manifest
    path "normalization_factors.tsv"

    // Log
    path "normalize_coverage_tracks.log", emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  # Ensure ERR trap propagates into functions/subshells (Bash)
  set -o errtrace
  export LC_ALL=C

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a normalize_coverage_tracks.log)
  exec 2> >(tee -a normalize_coverage_tracks.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "normalize_coverage_tracks" "Unexpected process failure" "Check normalize_coverage_tracks.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "NORMALIZE | START | sample=${sample_id} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="${sample_id}"
  CONDITION="${condition}"
  TIMEPOINT="${timepoint}"
  REPLICATE="${replicate}"
  THREADS=${task.cpus}
  
  COUNTS_MASTER="${counts_master_tsv}"
  GENOME_FA="${genome_fa}"
  GENES_TSV="${genes_tsv}"

  # Gene-end normalization (optional alternative CPM denominator; see
  # params.norm.gene_end_method doc in nextflow.config)
  GENE_END_METHOD="${params.norm?.gene_end_method ?: 'none'}"
  GENE_END_WINDOW=${params.norm?.gene_end_window ?: 500}
  GENE_END_MIN_READS=${params.norm?.gene_end_min_reads ?: 10}
  
  # Input bedGraphs
  POS3="${pos3_bg}"
  NEG3="${neg3_bg}"
  POS5="${pos5_bg}"
  NEG5="${neg5_bg}"
  AM3P_POS="${am3p_pos_bg}"
  AM3P_NEG="${am3p_neg_bg}"
  AM5P_POS="${am5p_pos_bg}"
  AM5P_NEG="${am5p_neg_bg}"

  # Feature flags
  EMIT_BW=\$([[ "${params.norm?.emit_bw}" == "false" ]] && echo 0 || echo 1)
  EMIT_SICPM=\$([[ "${params.norm?.emit_sicpm}" == "false" ]] && echo 0 || echo 1)
  EMIT_ALLMAP=\$([[ "${params.norm?.emit_allmap}" == "false" ]] && echo 0 || echo 1)
  FORCE_SORT=\$([[ "${params.force_sort_bedgraph}" == "true" ]] && echo 1 || echo 0)
  
  CONTROL_LABEL="${params.control_label ?: 'CTRL'}"
  TIMEOUT_BW=${params.norm?.timeout_bw ?: 900}

  # Auto-detect 5' track generation
  EMIT_5P_SETTING="${params.norm?.emit_5p}"
  if [[ "\${EMIT_5P_SETTING}" == "true" ]]; then
    EMIT_5P=1
  elif [[ "\${EMIT_5P_SETTING}" == "false" ]]; then
    EMIT_5P=0
  else
    # Auto mode: enable only if 5' inputs are NON-EMPTY. main.nf always supplies
    # an existing empty sentinel (EMPTY_5P_*.bedgraph), so testing -e here always
    # passed and forced 5' normalization of empty inputs every run. Test -s so the
    # intended "skip when there is no 5' data" actually triggers.
    if [[ -s "\${POS5}" || -s "\${NEG5}" ]]; then
      EMIT_5P=1
    else
      EMIT_5P=0
    fi
  fi

  echo "NORMALIZE | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "NORMALIZE | CONFIG | Condition: \${CONDITION}"
  echo "NORMALIZE | CONFIG | Timepoint: \${TIMEPOINT}"
  echo "NORMALIZE | CONFIG | Replicate: \${REPLICATE}"
  echo "NORMALIZE | CONFIG | Threads: \${THREADS}"
  echo "NORMALIZE | CONFIG | Control label: \${CONTROL_LABEL}"
  echo "NORMALIZE | CONFIG | Emit BigWig: \$([ \${EMIT_BW} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | Emit siCPM: \$([ \${EMIT_SICPM} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | Emit allMap: \$([ \${EMIT_ALLMAP} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | Emit 5' tracks: \$([ \${EMIT_5P} -eq 1 ] && echo "yes" || echo "no (auto)")"
  echo "NORMALIZE | CONFIG | Force bedGraph sort: \$([ \${FORCE_SORT} -eq 1 ] && echo "yes" || echo "no")"
  echo "NORMALIZE | CONFIG | BigWig timeout: \${TIMEOUT_BW}s"

  # Create output directories
  mkdir -p 3p 5p
  echo "NORMALIZE | CONFIG | Output directories created"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "NORMALIZE | VALIDATE | Checking input files..."

  # Validate counts master file
  if [[ ! -s "\${COUNTS_MASTER}" ]]; then
    tracktx_error "normalize_coverage_tracks" "Counts master file missing or empty: \${COUNTS_MASTER}" "Check quantify_reads_per_gene produced counts TSV"
  fi

  COUNTS_SIZE=\$(tracktx_size "\${COUNTS_MASTER}")
  COUNTS_LINES=\$(wc -l < "\${COUNTS_MASTER}" | tr -d ' ')
  echo "NORMALIZE | VALIDATE | Counts master: \${COUNTS_SIZE} bytes (\${COUNTS_LINES} lines)"

  # Validate input bedGraphs
  INPUT_COUNT=0
  for BG in "\${POS3}" "\${NEG3}" "\${AM3P_POS}" "\${AM3P_NEG}"; do
    if [[ -s "\${BG}" ]]; then
      INPUT_COUNT=\$((INPUT_COUNT + 1))
      BG_SIZE=\$(tracktx_size "\${BG}")
      echo "NORMALIZE | VALIDATE | Input bedGraph: \$(basename \${BG}) (\${BG_SIZE} bytes)"
    fi
  done

  echo "NORMALIZE | VALIDATE | Found \${INPUT_COUNT} non-empty input bedGraphs"

  # Shared resolver (bin/tracktx_error_fragment.sh): micromamba (container) ->
  # /opt/conda (container fallback) -> bare python3 (conda profile/local)
  tracktx_resolve_python

  # Validate tools
  if ! \${PYTHON_CMD} --version >/dev/null 2>&1; then
    tracktx_error "normalize_coverage_tracks" "Python not found (tried: \${PYTHON_CMD})" "Use -profile docker or install Python"
  fi
  echo "NORMALIZE | VALIDATE | python: \$(\${PYTHON_CMD} --version 2>&1)"
  if ! command -v awk >/dev/null 2>&1; then
    tracktx_error "normalize_coverage_tracks" "Required tool not found: awk" "Use -profile docker"
  fi
  echo "NORMALIZE | VALIDATE | awk: \$(which awk)"
  if [[ \${EMIT_BW} -eq 1 ]] && ! command -v bedGraphToBigWig >/dev/null 2>&1; then
    tracktx_error "normalize_coverage_tracks" "bedGraphToBigWig not found (required for BigWig)" "Install UCSC tools or use -profile docker"
  fi
  if [[ \${EMIT_BW} -eq 1 ]]; then
    echo "NORMALIZE | VALIDATE | bedGraphToBigWig: \$(which bedGraphToBigWig)"
  fi

  ###########################################################################
  # 3) COMPUTE NORMALIZATION FACTORS
  ###########################################################################

  echo "NORMALIZE | FACTORS | Computing CPM and siCPM scaling factors..."

  # Python script to compute factors from counts master
  \${PYTHON_CMD} - "\${COUNTS_MASTER}" "\${SAMPLE_ID}" "\${CONTROL_LABEL}" > factors.tmp <<'PYSCRIPT'
import sys
import csv

counts_file = sys.argv[1]
sample_id = sys.argv[2]
control_label = sys.argv[3].strip().lower()

# Read counts file
with open(counts_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    rows = list(reader)

# Helper function for case-insensitive comparison
def normalize_str(s):
    return (s or '').strip().lower()

# Find current sample
sample_row = None
for row in rows:
    if row.get('sample') == sample_id:
        sample_row = row
        break

if sample_row is None:
    print("0.0000000000\\t0.0000000000")
    sys.exit(0)

# Extract sample counts
sample_main = int(sample_row.get('main_reads', 0))
sample_spike = int(sample_row.get('spike_reads', 0))

# Compute CPM factor
if sample_main > 0:
    fac_cpm = 1_000_000.0 / sample_main
else:
    fac_cpm = 0.0

# Find control sample for siCPM.
#
# Match on condition == control_label and pick the LOWEST replicate number.
# Merged replicates are emitted with replicate = 0 (not 1), so treating
# replicate as a number and taking the minimum selects rep 1 when present
# and the merged rep-0 track when replicates are pooled. Only consider rows
# that actually have spike reads.
#
# IMPORTANT: there is intentionally NO fallback to "first sample with
# spike_reads > 0" when control_label doesn't match anything. An earlier
# version of this script did that silently, which — for any dataset using
# descriptive condition names instead of the literal default "CTRL" (i.e.
# every real dataset, unless control_label is set) — picked an arbitrary
# sample (often the treated condition, not the baseline) as the siCPM
# reference and produced a real-looking but wrong number with no warning.
control_row = None

def _rep_num(r):
    s = str(r.get('replicate', '')).strip().lower().lstrip('r')
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return 10**9

def _timepoint(r):
    try:
        return float(r.get('timepoint', ''))
    except (ValueError, TypeError):
        return None

control_candidates = [
    r for r in rows
    if normalize_str(r.get('condition', '')) == control_label
    and int(r.get('spike_reads', 0) or 0) > 0
]

# AUTO-DETECT fallback (replaces requiring control_label to be set by hand
# for the common case): if control_label didn't match anything, every
# TrackTx samplesheet already carries a structured 'timepoint' column where
# the untreated/baseline sample is the minimum timepoint in the timecourse
# (0 in every dataset used so far) -- independent of whatever free-text
# condition name a given study happens to use. This is trusted ONLY when
# it is unambiguous: exactly one distinct condition value present at that
# minimum timepoint with spike reads. If more than one condition shares
# the minimum timepoint (e.g. a "primed" and an "unprimed" arm both
# starting at timepoint 0), auto-detection refuses rather than guessing
# which arm this sample belongs to -- same never-guess-silently principle
# as the control_label match above, just narrowed to the genuinely
# ambiguous case instead of firing for every dataset that simply doesn't
# use the literal default control_label.
auto_detail = ""
if not control_candidates:
    tp_rows = [(r, _timepoint(r)) for r in rows]
    tp_rows = [(r, tp) for r, tp in tp_rows
               if tp is not None and int(r.get('spike_reads', 0) or 0) > 0]
    if tp_rows:
        min_tp = min(tp for _r, tp in tp_rows)
        at_min = [r for r, tp in tp_rows if tp == min_tp]
        distinct_conditions = sorted(set(normalize_str(r.get('condition', '')) for r in at_min))
        if len(distinct_conditions) == 1:
            control_candidates = at_min
            auto_detail = f"auto:timepoint={min_tp:g},condition={distinct_conditions[0]}"
        else:
            auto_detail = f"ambiguous:timepoint={min_tp:g},conditions={'|'.join(distinct_conditions)}"

used_auto = bool(auto_detail) and auto_detail.startswith("auto:")
if control_candidates:
    control_row = sorted(control_candidates, key=_rep_num)[0]

# Compute siCPM factor
if control_row and sample_spike > 0:
    control_main = int(control_row.get('main_reads', 0))
    control_spike = int(control_row.get('spike_reads', 0))

    if control_spike > 0 and control_main > 0:
        fac_sicpm = (control_spike / float(sample_spike)) * (1_000_000.0 / control_main)
        status = "ok_auto:" + auto_detail.split("auto:", 1)[1] if used_auto else "ok"
    else:
        fac_sicpm = 0.0
        status = "control_has_no_spike"
elif sample_spike <= 0:
    fac_sicpm = 0.0
    status = "sample_has_no_spike"
elif auto_detail:
    fac_sicpm = 0.0
    status = "no_control_label_match:" + auto_detail
else:
    fac_sicpm = 0.0
    status = "no_control_label_match"

print(f"{fac_cpm:.10f}\\t{fac_sicpm:.10f}\\t{status}")
PYSCRIPT

  # Read computed factors
  read -r FAC_CPM FAC_SICPM CONTROL_STATUS < factors.tmp

  echo "NORMALIZE | FACTORS | CPM factor: \${FAC_CPM}"
  echo "NORMALIZE | FACTORS | siCPM factor: \${FAC_SICPM}"

  # Validate CPM factor
  if awk -v x="\${FAC_CPM}" 'BEGIN{exit (x>0?0:1)}'; then
    echo "NORMALIZE | FACTORS | CPM normalization enabled"
  else
    tracktx_error "normalize_coverage_tracks" "Cannot compute CPM (sample reads = 0)" "Check quantify_reads_per_gene output and counts TSV"
  fi

  # Check siCPM availability
  if awk -v x="\${FAC_SICPM}" 'BEGIN{exit (x>0?0:1)}'; then
    case "\${CONTROL_STATUS}" in
      ok_auto:*)
        echo "NORMALIZE | FACTORS | siCPM normalization enabled (control auto-detected: \${CONTROL_STATUS#ok_auto:})"
        ;;
      *)
        echo "NORMALIZE | FACTORS | siCPM normalization enabled"
        ;;
    esac
    SICPM_AVAILABLE=1
  else
    echo "NORMALIZE | FACTORS | WARNING: siCPM disabled (factor = 0)"
    case "\${CONTROL_STATUS}" in
      no_control_label_match:ambiguous:*)
        echo "NORMALIZE | FACTORS | Cause: no sample has condition == [\${CONTROL_LABEL}] (params.control_label), and automatic timepoint-based detection found MORE THAN ONE condition at the minimum timepoint, so it refused to guess: \${CONTROL_STATUS#no_control_label_match:ambiguous:}"
        echo "NORMALIZE | FACTORS | Fix: set control_label in your params file to the EXACT baseline condition name for the arm this sample belongs to."
        ;;
      no_control_label_match)
        echo "NORMALIZE | FACTORS | Cause: no sample has condition == [\${CONTROL_LABEL}] (params.control_label) with spike_reads > 0, and automatic timepoint-based detection found no usable timepoint data either."
        echo "NORMALIZE | FACTORS | Fix: set control_label in your params file to the EXACT baseline condition name used in your samplesheet (e.g. control_label: no_heat_shock), or ensure the samplesheet's timepoint column is populated."
        ;;
      sample_has_no_spike)
        echo "NORMALIZE | FACTORS | Cause: this sample (\${SAMPLE_ID}) itself has 0 spike-in reads."
        ;;
      control_has_no_spike)
        echo "NORMALIZE | FACTORS | Cause: the matched control sample has 0 spike-in or main reads."
        ;;
      *)
        echo "NORMALIZE | FACTORS | Possible causes: no spike-in reads, no control sample"
        ;;
    esac
    SICPM_AVAILABLE=0
  fi

  ###########################################################################
  # 3b) COMPUTE GENE-END NORMALIZATION FACTOR (optional)
  ###########################################################################
  #
  # Alternative CPM denominator: total raw 3' signal inside a per-gene
  # gene-end region (tes_window: last GENE_END_WINDOW bp of the gene body
  # ending at the TES, on the gene's own strand; gene_body: the whole gene
  # footprint), summed only over genes clearing GENE_END_MIN_READS. Written
  # to normalization_factors.tsv as gene_end_<method> -- not applied to any
  # bedGraph/BigWig track (CPM/siCPM remain the only scaled track outputs).
  # Uses POS3/NEG3 (raw, pre-scaling 3' bedGraphs), which module 06 already
  # writes coordinate-sorted, matching the sortedness assumption the rest of
  # this script makes for BigWig conversion.

  FAC_GENEEND="0.0000000000"
  if [[ "\${GENE_END_METHOD}" != "none" ]]; then
    echo "NORMALIZE | FACTORS | Computing gene-end factor (method=\${GENE_END_METHOD}, window=\${GENE_END_WINDOW}bp, min_reads=\${GENE_END_MIN_READS})..."

    if [[ ! -s "\${GENES_TSV}" ]]; then
      echo "NORMALIZE | FACTORS | WARNING: gene-end factor disabled (genes.tsv missing/empty)"
    else
      # genes.tsv: gene_id  gene_name  chr  strand  start  end  tss  tes  biotype
      # (start/end/tss/tes are 1-based GTF coords; emit 0-based BED6 here)
      awk -F'\t' -v OFS='\t' -v method="\${GENE_END_METHOD}" -v win="\${GENE_END_WINDOW}" '
        NR==1 { next }
        {
          chrom=\$3; strand=\$4; start=\$5; end=\$6; tes=\$8
          if (method == "gene_body") {
            lo = start - 1; hi = end
          } else if (strand == "+") {
            lo = tes - win; if (lo < 0) lo = 0; hi = tes
          } else {
            lo = tes - 1; hi = (tes - 1) + win
          }
          if (hi <= lo) next
          print chrom, lo, hi, \$1, 0, strand
        }
      ' "\${GENES_TSV}" | LC_ALL=C sort -k1,1 -k2,2n > gene_end_regions.bed

      GENE_END_TOTAL=0
      if [[ -s gene_end_regions.bed ]]; then
        awk -F'\t' '\$6=="+"' gene_end_regions.bed > gene_end_regions.pos.bed
        awk -F'\t' '\$6=="-"' gene_end_regions.bed > gene_end_regions.neg.bed

        if [[ -s gene_end_regions.pos.bed && -s "\${POS3}" ]]; then
          POS_SUM=\$(bedtools map -a gene_end_regions.pos.bed -b "\${POS3}" -c 4 -o sum -null 0 \\
            | awk -v minr="\${GENE_END_MIN_READS}" '{v=\$7; if (v<0) v=-v; if (v>=minr) s+=v} END{printf "%.0f", s+0}')
          GENE_END_TOTAL=\$(awk -v a="\${GENE_END_TOTAL}" -v b="\${POS_SUM:-0}" 'BEGIN{printf "%.0f", a+b}')
        fi
        if [[ -s gene_end_regions.neg.bed && -s "\${NEG3}" ]]; then
          NEG_SUM=\$(bedtools map -a gene_end_regions.neg.bed -b "\${NEG3}" -c 4 -o sum -null 0 \\
            | awk -v minr="\${GENE_END_MIN_READS}" '{v=\$7; if (v<0) v=-v; if (v>=minr) s+=v} END{printf "%.0f", s+0}')
          GENE_END_TOTAL=\$(awk -v a="\${GENE_END_TOTAL}" -v b="\${NEG_SUM:-0}" 'BEGIN{printf "%.0f", a+b}')
        fi
        rm -f gene_end_regions.pos.bed gene_end_regions.neg.bed
      fi
      rm -f gene_end_regions.bed

      if awk -v x="\${GENE_END_TOTAL}" 'BEGIN{exit (x>0?0:1)}'; then
        FAC_GENEEND=\$(awk -v t="\${GENE_END_TOTAL}" 'BEGIN{printf "%.10f", 1000000.0/t}')
        echo "NORMALIZE | FACTORS | Gene-end factor: \${FAC_GENEEND} (total=\${GENE_END_TOTAL} reads across qualifying genes)"
      else
        echo "NORMALIZE | FACTORS | WARNING: gene-end factor disabled (no gene cleared gene_end_min_reads=\${GENE_END_MIN_READS})"
      fi
    fi
  fi

  ###########################################################################
  # 4) PREPARE GENOME SIZES
  ###########################################################################

  echo "NORMALIZE | GENOME | Preparing chromosome sizes..."

  # Check for existing .fai
  GENOME_FA_SRC="${genome_fa.toString()}"
  if [[ -s "\${GENOME_FA_SRC}.fai" && ! -s "\${GENOME_FA}.fai" ]]; then
    ln -sf "\${GENOME_FA_SRC}.fai" "\${GENOME_FA}.fai" 2>/dev/null || true
  fi

  # Create genome.sizes
  if [[ -s "\${GENOME_FA}.fai" ]]; then
    echo "NORMALIZE | GENOME | Using existing FASTA index"
    cut -f1,2 "\${GENOME_FA}.fai" > genome.sizes
  else
    echo "NORMALIZE | GENOME | Creating FASTA index..."
    samtools faidx "\${GENOME_FA}"
    cut -f1,2 "\${GENOME_FA}.fai" > genome.sizes
  fi

  # Sort and validate genome.sizes
  LC_ALL=C sort -k1,1 -u -o genome.sizes genome.sizes

  CHR_COUNT=\$(wc -l < genome.sizes | tr -d ' ')
  TOTAL_SIZE=\$(awk '{sum+=\$2} END{print sum}' genome.sizes)
  echo "NORMALIZE | GENOME | Chromosomes: \${CHR_COUNT}"
  echo "NORMALIZE | GENOME | Total size: \${TOTAL_SIZE} bp"

  ###########################################################################
  # 5) HELPER FUNCTIONS
  ###########################################################################

  # Convert bedGraph to BigWig
  make_bigwig() {
    local bedgraph="\$1"
    local bigwig="\$2"
    
    if [[ ! -s "\${bedgraph}" ]]; then
      echo "NORMALIZE | BIGWIG | Empty input, creating empty BigWig: \$(basename \${bigwig})"
      : > "\${bigwig}"
      return 0
    fi
    
    if [[ \${EMIT_BW} -eq 0 ]]; then
      : > "\${bigwig}"
      return 0
    fi
    
    echo "NORMALIZE | BIGWIG | Converting: \$(basename \${bedgraph}) → \$(basename \${bigwig})"
    
    # Optional sorting
    if [[ \${FORCE_SORT} -eq 1 ]]; then
      echo "NORMALIZE | BIGWIG | Sorting bedGraph..."
      # Give sort most of the task's RAM so it stays in memory instead of spilling
      # thousands of tiny temp files to disk (catastrophic on slow/USB work dirs).
      # Falls back to disk only if truly needed, using a fast temp dir (never the
      # USB-backed work dir via "-T ."). Override with SORT_MEM / SORT_TMPDIR.
      : "\${SORT_MEM:=\$(( ${task.memory.toGiga()} * 70 / 100 ))G}"
      : "\${SORT_TMP:=\${SORT_TMPDIR:-/tmp}}"
      mkdir -p "\${SORT_TMP}" 2>/dev/null || SORT_TMP=/tmp
      if ! LC_ALL=C sort -S "\${SORT_MEM}" -T "\${SORT_TMP}" --parallel="\${SORT_PARALLEL:-\${THREADS}}" -k1,1 -k2,2n "\${bedgraph}" > "\${bedgraph}.sorted"; then
        echo "NORMALIZE | ERROR | sort failed (likely OOM) for: \${bedgraph}"
        rm -f "\${bedgraph}.sorted"
        return 1
      fi
      mv -f "\${bedgraph}.sorted" "\${bedgraph}"
    fi
    
    # Count lines
    LINE_COUNT=\$(awk 'BEGIN{n=0} \$0!~/^(track|browser|#)/{n++} END{print n}' "\${bedgraph}")
    echo "NORMALIZE | BIGWIG | bedGraph lines: \${LINE_COUNT}"
    
    # Convert with timeout
    if timeout "\${TIMEOUT_BW}" bedGraphToBigWig "\${bedgraph}" genome.sizes "\${bigwig}"; then
      BW_SIZE=\$(tracktx_size "\${bigwig}")
      echo "NORMALIZE | BIGWIG | Created: \$(basename \${bigwig}) (\${BW_SIZE} bytes)"
    else
      echo "NORMALIZE | BIGWIG | WARNING: Conversion failed or timed out, creating empty BigWig"
      : > "\${bigwig}"
    fi
  }

  # Normalize bedGraph (single-pass CPM + siCPM)
  normalize_bedgraph() {
    local input_bg="\$1"
    local end_label="\$2"      # "3p" or "5p"
    local strand="\$3"          # "pos" or "neg"
    local set_label="\$4"       # "main" or "allMap"
    
    # Build output file names
    local set_suffix=""
    if [[ "\${set_label}" == "allMap" ]]; then
      set_suffix=".allMap"
    fi
    
    local prefix="\${SAMPLE_ID}\${set_suffix}.\${end_label}.\${strand}"
    local out_cpm_bg="\${end_label}/\${prefix}.cpm.bedgraph"
    local out_cpm_bw="\${end_label}/\${prefix}.cpm.bw"
    local out_sicpm_bg="\${end_label}/\${prefix}.sicpm.bedgraph"
    local out_sicpm_bw="\${end_label}/\${prefix}.sicpm.bw"
    
    echo "NORMALIZE | SCALE | Processing: \${end_label} \${set_label} \${strand}"
    
    if [[ ! -s "\${input_bg}" ]]; then
      echo "NORMALIZE | SCALE | WARNING: Input empty, creating empty outputs"
      : > "\${out_cpm_bg}"
      : > "\${out_cpm_bw}"
      : > "\${out_sicpm_bg}"
      : > "\${out_sicpm_bw}"
      return 0
    fi
    
    # Count input lines
    INPUT_LINES=\$(awk 'BEGIN{n=0} \$0!~/^(track|browser|#)/{n++} END{print n}' "\${input_bg}")
    echo "NORMALIZE | SCALE | Input lines: \${INPUT_LINES}"
    
    # Normalize with awk (single pass for both CPM and siCPM)
    START_TIME=\$(date +%s)
    
    if [[ \${EMIT_SICPM} -eq 1 && \${SICPM_AVAILABLE} -eq 1 ]]; then
      echo "NORMALIZE | SCALE | Writing CPM and siCPM..."
      awk -v fc="\${FAC_CPM}" -v fs="\${FAC_SICPM}" -v OFS='\\t' '
        BEGIN { OFMT = "%.10f" }
        (NF>=4) && (\$0!~/^(track|browser|#)/) {
          cpm_val = \$4 * fc
          sicpm_val = \$4 * fs
          print \$1, \$2, \$3, cpm_val > "'"\${out_cpm_bg}"'"
          print \$1, \$2, \$3, sicpm_val > "'"\${out_sicpm_bg}"'"
        }
      ' "\${input_bg}"
    else
      echo "NORMALIZE | SCALE | Writing CPM only..."
      awk -v fc="\${FAC_CPM}" -v OFS='\\t' '
        BEGIN { OFMT = "%.10f" }
        (NF>=4) && (\$0!~/^(track|browser|#)/) {
          print \$1, \$2, \$3, \$4 * fc
        }
      ' "\${input_bg}" > "\${out_cpm_bg}"
      : > "\${out_sicpm_bg}"
    fi
    
    END_TIME=\$(date +%s)
    ELAPSED=\$((END_TIME - START_TIME))
    if [[ \${ELAPSED} -lt 1 ]]; then ELAPSED=1; fi
    RATE=\$((INPUT_LINES / ELAPSED))
    
    echo "NORMALIZE | SCALE | Normalization complete in \${ELAPSED}s (~\${RATE} lines/s)"
    
    # Report output sizes
    CPM_LINES=\$(wc -l < "\${out_cpm_bg}" 2>/dev/null | tr -d ' ' || echo 0)
    CPM_SIZE=\$(tracktx_size "\${out_cpm_bg}")
    echo "NORMALIZE | SCALE | CPM bedGraph: \${CPM_LINES} lines (\${CPM_SIZE} bytes)"
    
    if [[ -s "\${out_sicpm_bg}" ]]; then
      SICPM_LINES=\$(wc -l < "\${out_sicpm_bg}" 2>/dev/null | tr -d ' ' || echo 0)
      SICPM_SIZE=\$(tracktx_size "\${out_sicpm_bg}")
      echo "NORMALIZE | SCALE | siCPM bedGraph: \${SICPM_LINES} lines (\${SICPM_SIZE} bytes)"
    fi
    
    # Convert to BigWig
    make_bigwig "\${out_cpm_bg}" "\${out_cpm_bw}"
    if [[ -s "\${out_sicpm_bg}" ]]; then
      make_bigwig "\${out_sicpm_bg}" "\${out_sicpm_bw}"
    else
      : > "\${out_sicpm_bw}"
    fi
    
    # Add to manifest. Written to a per-call fragment file (not appended to
    # tracks_manifest.tsv directly) since normalize_bedgraph now runs several
    # at once in the parallel fan-out below -- concurrent appends to one file
    # from background jobs aren't safe. Fragments are concatenated after all
    # jobs finish.
    local manifest_frag="manifest_\${set_label}_\${end_label}_\${strand}.tsv"
    : > "\${manifest_frag}"
    echo -e "\${SAMPLE_ID}\\t\${end_label}\\t\${set_label}\\t\${strand}\\tcpm\\t\${out_cpm_bg}" >> "\${manifest_frag}"
    if [[ -s "\${out_sicpm_bg}" ]]; then
      echo -e "\${SAMPLE_ID}\\t\${end_label}\\t\${set_label}\\t\${strand}\\tsicpm\\t\${out_sicpm_bg}" >> "\${manifest_frag}"
    fi
    
    echo "NORMALIZE | SCALE | Complete: \${end_label} \${set_label} \${strand}"
  }

  ###########################################################################
  # 6-9) NORMALIZE ALL ENABLED TRACK SETS (parallel fan-out)
  ###########################################################################
  # Up to 8 independent normalize_bedgraph calls (main/allMap × 3p/5p ×
  # pos/neg); each writes distinct output files, so they're safe to run
  # concurrently. Throttled the same way module 06 throttles its coverage
  # jobs: cap concurrency by task memory (each job's sort budgeted at ~2GB
  # peak) and never more than the allotted CPUs, then divide SORT_MEM/THREADS
  # across however many run at once so concurrent sorts don't oversubscribe
  # either.
  echo "NORMALIZE | SCALE | Normalizing enabled track sets (parallel)..."

  # Initialize manifest
  : > tracks_manifest.tsv

  MEM_GB=${task.memory.toGiga()}
  MAX_PAR=\$(( MEM_GB * 70 / 100 / 2 ))
  [ "\${MAX_PAR}" -lt 1 ] && MAX_PAR=1
  [ "\${MAX_PAR}" -gt "\${THREADS}" ] && MAX_PAR=\${THREADS}
  MAX_PAR=\${NORMALIZE_MAX_PAR:-\${MAX_PAR}}

  SORT_MEM_GB=\$(( MEM_GB * 50 / 100 / MAX_PAR ))
  [ "\${SORT_MEM_GB}" -lt 1 ] && SORT_MEM_GB=1
  export SORT_MEM="\${SORT_MEM_GB}G"
  export SORT_TMP="\${SORT_TMPDIR:-/tmp}"
  SORT_PARALLEL=\$(( THREADS / MAX_PAR ))
  [ "\${SORT_PARALLEL}" -lt 1 ] && SORT_PARALLEL=1
  export SORT_PARALLEL
  echo "NORMALIZE | SCALE | Concurrency: \${MAX_PAR} parallel job(s), SORT_MEM=\${SORT_MEM} each, SORT_PARALLEL=\${SORT_PARALLEL} (task mem=\${MEM_GB}G)"

  NORM_FAIL_FLAG="normalize_fail.flag"
  rm -f "\${NORM_FAIL_FLAG}"

  NORM_PIDS=()
  launch_norm() {
    # Block until fewer than MAX_PAR of OUR jobs are still alive. Same
    # tracked-PID pattern as module 06's launch_cov -- a bare \`wait\` would
    # hang on the \`tee\` process-substitution children from \`exec >\` above.
    while :; do
      local alive=0 p
      for p in "\${NORM_PIDS[@]:-}"; do
        [ -n "\${p}" ] && kill -0 "\${p}" 2>/dev/null && alive=\$(( alive + 1 ))
      done
      [ "\${alive}" -lt "\${MAX_PAR}" ] && break
      sleep 0.5
    done
    ( normalize_bedgraph "\$1" "\$2" "\$3" "\$4" || echo "FAIL: \$2 \$3 \$4" >> "\${NORM_FAIL_FLAG}" ) &
    NORM_PIDS+=(\$!)
  }

  launch_norm "\${POS3}" "3p" "pos" "main"
  launch_norm "\${NEG3}" "3p" "neg" "main"

  if [[ \${EMIT_5P} -eq 1 ]]; then
    launch_norm "\${POS5}" "5p" "pos" "main"
    launch_norm "\${NEG5}" "5p" "neg" "main"
  else
    echo "NORMALIZE | MAIN5P | Skipping 5' tracks (not enabled)"
  fi

  if [[ \${EMIT_ALLMAP} -eq 1 ]]; then
    launch_norm "\${AM3P_POS}" "3p" "pos" "allMap"
    launch_norm "\${AM3P_NEG}" "3p" "neg" "allMap"
  else
    echo "NORMALIZE | ALLMAP3P | Skipping allMap 3' tracks (not enabled)"
  fi

  if [[ \${EMIT_ALLMAP} -eq 1 && \${EMIT_5P} -eq 1 ]]; then
    launch_norm "\${AM5P_POS}" "5p" "pos" "allMap"
    launch_norm "\${AM5P_NEG}" "5p" "neg" "allMap"
  else
    echo "NORMALIZE | ALLMAP5P | Skipping allMap 5' tracks"
  fi

  for p in "\${NORM_PIDS[@]}"; do
    wait "\${p}" || true
  done

  if [[ -s "\${NORM_FAIL_FLAG}" ]]; then
    echo "NORMALIZE | ERROR | Failed normalization jobs:"
    sed 's/^/NORMALIZE | ERROR |   /' "\${NORM_FAIL_FLAG}"
    tracktx_error "normalize_coverage_tracks" "One or more track normalization jobs failed" "Check normalize_coverage_tracks.log for per-job error messages"
  fi

  # Merge per-job manifest fragments (see normalize_bedgraph) into the real
  # manifest, then clean them up.
  cat manifest_*.tsv >> tracks_manifest.tsv 2>/dev/null || true
  rm -f manifest_*.tsv

  echo "NORMALIZE | SCALE | All enabled track sets complete"

  ###########################################################################
  # 9c) GUARANTEE PER-SAMPLE-REPORT BIGWIGS EXIST
  ###########################################################################
  # The per-sample report (module 14) receives these 4 CPM BigWigs through the
  # report_bw channel, so its track-link availability is decided from real
  # produced files instead of probing the publish directory. Touch a 0-byte
  # placeholder for any not produced (allMap when emit_allmap=false; all .bw when
  # emit_bw=false) so the output binding always succeeds; 0-byte = not available.
  for RBW in \\
    "3p/\${SAMPLE_ID}.3p.pos.cpm.bw" \\
    "3p/\${SAMPLE_ID}.3p.neg.cpm.bw" \\
    "3p/\${SAMPLE_ID}.allMap.3p.pos.cpm.bw" \\
    "3p/\${SAMPLE_ID}.allMap.3p.neg.cpm.bw"; do
    [[ -e "\${RBW}" ]] || : > "\${RBW}"
  done

  ###########################################################################
  # 9d) GUARANTEE MAIN 5' CPM/siCPM BEDGRAPHS EXIST (cpm5p_bg/sicpm5p_bg outputs)
  ###########################################################################
  # GRO-seq routes its functional-region and Pol-II density signal through the
  # 5' end (params.signal_end=5p; see main.nf STEP 11/12), so these need to be
  # real declared outputs, not just files that happen to land in 5p/ when
  # EMIT_5P=1. When EMIT_5P=0 (PRO-seq default, or emit_5p=false) section 7
  # above never ran and these paths don't exist yet -- touch empty
  # placeholders so the output binding always succeeds; an empty bedgraph is
  # already the established "not available" convention (module 06's sentinel
  # EMPTY_5P_*.bedgraph, EMIT_SICPM=0 path above, etc.).
  for MAIN5P in \\
    "5p/\${SAMPLE_ID}.5p.pos.cpm.bedgraph" \\
    "5p/\${SAMPLE_ID}.5p.neg.cpm.bedgraph" \\
    "5p/\${SAMPLE_ID}.5p.pos.sicpm.bedgraph" \\
    "5p/\${SAMPLE_ID}.5p.neg.sicpm.bedgraph"; do
    [[ -e "\${MAIN5P}" ]] || : > "\${MAIN5P}"
  done

  ###########################################################################
  # 10) CREATE NORMALIZATION FACTORS FILE
  ###########################################################################

  echo "NORMALIZE | OUTPUT | Writing normalization factors..."

  cat > normalization_factors.tsv <<FACTOREOF
method  factor
CPM \${FAC_CPM}
siCPM \${FAC_SICPM}
FACTOREOF
  if [[ "\${GENE_END_METHOD}" != "none" ]]; then
    echo -e "gene_end_\${GENE_END_METHOD}\t\${FAC_GENEEND}" >> normalization_factors.tsv
  fi

  echo "NORMALIZE | OUTPUT | Normalization factors written"

  ###########################################################################
  # 11) CREATE LEGACY SYMLINKS
  ###########################################################################

  echo "NORMALIZE | LEGACY | Creating legacy symlinks for compatibility..."

  # Legacy hard links for downstream compatibility (CIFS mount: no symlink support)
  ln -f "3p/\${SAMPLE_ID}.3p.pos.cpm.bedgraph" "3p/\${SAMPLE_ID}_pos3_cpm.bedgraph"
  ln -f "3p/\${SAMPLE_ID}.3p.neg.cpm.bedgraph" "3p/\${SAMPLE_ID}_neg3_cpm.bedgraph"

  echo "NORMALIZE | LEGACY | Legacy symlinks created"

  ###########################################################################
  # 12) CREATE README
  ###########################################################################

  echo "NORMALIZE | README | Creating documentation..."

  cat > README_normalization.txt <<DOCEOF
NORMALIZED TRACKS — ${sample_id}
────────────────────────────────────────────────────────────────────────────
  CPM   = (raw_count / sample_reads) x 1,000,000
  siCPM = (raw_count / sample_spike) x (control_spike / control_reads) x 1,000,000

  Control for siCPM: condition == "\${CONTROL_LABEL}" (params.control_label),
  lowest replicate number preferred. MUST be set per dataset to your real
  baseline condition -- no fallback to "first sample with spike reads"; a
  non-matching label disables siCPM (factor = 0) rather than guessing.
  This sample: CPM factor=\${FAC_CPM}  siCPM factor=\${FAC_SICPM}
  \$([ \${SICPM_AVAILABLE} -eq 0 ] && echo "  -> siCPM DISABLED for this sample" || echo "  -> siCPM available")

  3p/${sample_id}.3p.{pos,neg}.{cpm,sicpm}.{bedgraph,bw}   — always generated
  5p/${sample_id}.5p.{pos,neg}.{cpm,sicpm}.{bedgraph,bw}   — always generated
  3p/${sample_id}.allMap.3p.*  (+ 5p equivalents)          — if emit_allmap
  normalization_factors.tsv, tracks_manifest.tsv (sample/end/set/strand/scale/path)

  BigWig timeout \${TIMEOUT_BW}s; force-sort bedGraph: \$([ \${FORCE_SORT} -eq 1 ] && echo "yes" || echo "no").
  Downstream: divergent-TX calling and functional-region assignment use the
  MAIN (not allMap) tracks on whichever end params.signal_end selects (3' for
  PRO-seq, 5' for GRO-seq); Pol-II density uses siCPM when available, else CPM.
DOCEOF

  echo "NORMALIZE | README | Documentation created"

  ###########################################################################
  # 13) VALIDATION AND SUMMARY
  ###########################################################################

  echo "NORMALIZE | VALIDATE | Validating outputs..."

  # Count output files
  CPM_BG_COUNT=\$(find 3p 5p -name "*.cpm.bedgraph" -type f 2>/dev/null | wc -l | tr -d ' ')
  SICPM_BG_COUNT=\$(find 3p 5p -name "*.sicpm.bedgraph" -type f 2>/dev/null | wc -l | tr -d ' ')
  CPM_BW_COUNT=\$(find 3p 5p -name "*.cpm.bw" -type f 2>/dev/null | wc -l | tr -d ' ')
  SICPM_BW_COUNT=\$(find 3p 5p -name "*.sicpm.bw" -type f 2>/dev/null | wc -l | tr -d ' ')
  MANIFEST_LINES=\$(wc -l < tracks_manifest.tsv | tr -d ' ')

  echo "NORMALIZE | VALIDATE | CPM bedGraphs: \${CPM_BG_COUNT}"
  echo "NORMALIZE | VALIDATE | siCPM bedGraphs: \${SICPM_BG_COUNT}"
  echo "NORMALIZE | VALIDATE | CPM BigWigs: \${CPM_BW_COUNT}"
  echo "NORMALIZE | VALIDATE | siCPM BigWigs: \${SICPM_BW_COUNT}"
  echo "NORMALIZE | VALIDATE | Manifest entries: \${MANIFEST_LINES}"

  # Check critical files exist
  for file in \\
    "3p/\${SAMPLE_ID}.3p.pos.cpm.bedgraph" \\
    "3p/\${SAMPLE_ID}.3p.neg.cpm.bedgraph" \\
    "normalization_factors.tsv" \\
    "tracks_manifest.tsv"; do
    
    if [[ ! -s "\${file}" ]]; then
      tracktx_error "normalize_coverage_tracks" "Missing or empty critical file: \${file}" "Check normalize_coverage_tracks.log in work dir"
    fi
  done


  echo "NORMALIZE | VALIDATE | All critical files present"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  # Calculate total output size
  TOTAL_SIZE=\$(du -sh . 2>/dev/null | cut -f1 || echo "unknown")

  echo "────────────────────────────────────────────────────────────────────────"
  echo "NORMALIZE | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "NORMALIZE | SUMMARY | CPM factor: \${FAC_CPM}"
  echo "NORMALIZE | SUMMARY | siCPM factor: \${FAC_SICPM}"
  echo "NORMALIZE | SUMMARY | CPM bedGraphs: \${CPM_BG_COUNT}"
  echo "NORMALIZE | SUMMARY | siCPM bedGraphs: \${SICPM_BG_COUNT}"
  echo "NORMALIZE | SUMMARY | BigWig files: \$((CPM_BW_COUNT + SICPM_BW_COUNT))"
  echo "NORMALIZE | SUMMARY | Total output: \${TOTAL_SIZE}"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "NORMALIZE | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}