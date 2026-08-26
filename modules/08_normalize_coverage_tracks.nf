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
//         condition, timepoint, replicate, counts_master_tsv, genes_unused)
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
          val(genes_unused)
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

    // Full directory trees
    path "3p/**", emit: tree3p
    path "5p/**", optional: true, emit: tree5p

    // Log
    path "normalize_coverage_tracks.log", emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a normalize_coverage_tracks.log)
  exec 2> >(tee -a normalize_coverage_tracks.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'tracktx_error "normalize_coverage_tracks" "Unexpected process failure" "Check normalize_coverage_tracks.log in work dir"' ERR

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

  COUNTS_SIZE=\$(stat -c%s "\${COUNTS_MASTER}" 2>/dev/null || stat -f%z "\${COUNTS_MASTER}" 2>/dev/null || echo "unknown")
  COUNTS_LINES=\$(wc -l < "\${COUNTS_MASTER}" | tr -d ' ')
  echo "NORMALIZE | VALIDATE | Counts master: \${COUNTS_SIZE} bytes (\${COUNTS_LINES} lines)"

  # Validate input bedGraphs
  INPUT_COUNT=0
  for BG in "\${POS3}" "\${NEG3}" "\${AM3P_POS}" "\${AM3P_NEG}"; do
    if [[ -s "\${BG}" ]]; then
      INPUT_COUNT=\$((INPUT_COUNT + 1))
      BG_SIZE=\$(stat -c%s "\${BG}" 2>/dev/null || stat -f%z "\${BG}" 2>/dev/null || echo "unknown")
      echo "NORMALIZE | VALIDATE | Input bedGraph: \$(basename \${BG}) (\${BG_SIZE} bytes)"
    fi
  done

  echo "NORMALIZE | VALIDATE | Found \${INPUT_COUNT} non-empty input bedGraphs"

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  else
    PYTHON_CMD="python3"
  fi

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
      if ! LC_ALL=C sort -S "\${SORT_MEM}" -T "\${SORT_TMP}" -k1,1 -k2,2n "\${bedgraph}" > "\${bedgraph}.sorted"; then
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
      BW_SIZE=\$(stat -c%s "\${bigwig}" 2>/dev/null || stat -f%z "\${bigwig}" 2>/dev/null || echo "unknown")
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
    CPM_SIZE=\$(stat -c%s "\${out_cpm_bg}" 2>/dev/null || stat -f%z "\${out_cpm_bg}" 2>/dev/null || echo "unknown")
    echo "NORMALIZE | SCALE | CPM bedGraph: \${CPM_LINES} lines (\${CPM_SIZE} bytes)"
    
    if [[ -s "\${out_sicpm_bg}" ]]; then
      SICPM_LINES=\$(wc -l < "\${out_sicpm_bg}" 2>/dev/null | tr -d ' ' || echo 0)
      SICPM_SIZE=\$(stat -c%s "\${out_sicpm_bg}" 2>/dev/null || stat -f%z "\${out_sicpm_bg}" 2>/dev/null || echo "unknown")
      echo "NORMALIZE | SCALE | siCPM bedGraph: \${SICPM_LINES} lines (\${SICPM_SIZE} bytes)"
    fi
    
    # Convert to BigWig
    make_bigwig "\${out_cpm_bg}" "\${out_cpm_bw}"
    if [[ -s "\${out_sicpm_bg}" ]]; then
      make_bigwig "\${out_sicpm_bg}" "\${out_sicpm_bw}"
    else
      : > "\${out_sicpm_bw}"
    fi
    
    # Add to manifest
    echo -e "\${SAMPLE_ID}\\t\${end_label}\\t\${set_label}\\t\${strand}\\tcpm\\t\${out_cpm_bg}" >> tracks_manifest.tsv
    if [[ -s "\${out_sicpm_bg}" ]]; then
      echo -e "\${SAMPLE_ID}\\t\${end_label}\\t\${set_label}\\t\${strand}\\tsicpm\\t\${out_sicpm_bg}" >> tracks_manifest.tsv
    fi
    
    echo "NORMALIZE | SCALE | Complete: \${end_label} \${set_label} \${strand}"
  }

  ###########################################################################
  # 6) NORMALIZE MAIN 3' TRACKS
  ###########################################################################

  echo "NORMALIZE | MAIN3P | Normalizing main 3' tracks..."

  # Initialize manifest
  : > tracks_manifest.tsv

  normalize_bedgraph "\${POS3}" "3p" "pos" "main"
  normalize_bedgraph "\${NEG3}" "3p" "neg" "main"

  echo "NORMALIZE | MAIN3P | Main 3' tracks complete"

  ###########################################################################
  # 7) NORMALIZE MAIN 5' TRACKS (if enabled)
  ###########################################################################

  if [[ \${EMIT_5P} -eq 1 ]]; then
    echo "NORMALIZE | MAIN5P | Normalizing main 5' tracks..."
    
    normalize_bedgraph "\${POS5}" "5p" "pos" "main"
    normalize_bedgraph "\${NEG5}" "5p" "neg" "main"
    
    echo "NORMALIZE | MAIN5P | Main 5' tracks complete"
  else
    echo "NORMALIZE | MAIN5P | Skipping 5' tracks (not enabled)"
  fi

  ###########################################################################
  # 8) NORMALIZE ALLMAP 3' TRACKS (if enabled)
  ###########################################################################

  if [[ \${EMIT_ALLMAP} -eq 1 ]]; then
    echo "NORMALIZE | ALLMAP3P | Normalizing allMap 3' tracks..."
    
    normalize_bedgraph "\${AM3P_POS}" "3p" "pos" "allMap"
    normalize_bedgraph "\${AM3P_NEG}" "3p" "neg" "allMap"
    
    echo "NORMALIZE | ALLMAP3P | AllMap 3' tracks complete"
  else
    echo "NORMALIZE | ALLMAP3P | Skipping allMap 3' tracks (not enabled)"
  fi

  ###########################################################################
  # 9) NORMALIZE ALLMAP 5' TRACKS (if both enabled)
  ###########################################################################

  if [[ \${EMIT_ALLMAP} -eq 1 && \${EMIT_5P} -eq 1 ]]; then
    echo "NORMALIZE | ALLMAP5P | Normalizing allMap 5' tracks..."
    
    normalize_bedgraph "\${AM5P_POS}" "5p" "pos" "allMap"
    normalize_bedgraph "\${AM5P_NEG}" "5p" "neg" "allMap"
    
    echo "NORMALIZE | ALLMAP5P | AllMap 5' tracks complete"
  else
    echo "NORMALIZE | ALLMAP5P | Skipping allMap 5' tracks"
  fi

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
  # 10) CREATE NORMALIZATION FACTORS FILE
  ###########################################################################

  echo "NORMALIZE | OUTPUT | Writing normalization factors..."

  cat > normalization_factors.tsv <<FACTOREOF
method  factor
CPM \${FAC_CPM}
siCPM \${FAC_SICPM}
FACTOREOF

  echo "NORMALIZE | OUTPUT | Normalization factors written"

  ###########################################################################
  # 11) CREATE LEGACY SYMLINKS
  ###########################################################################

  echo "NORMALIZE | LEGACY | Creating legacy symlinks for compatibility..."

  # Legacy symlinks for downstream compatibility
  ln -sf "\${SAMPLE_ID}.3p.pos.cpm.bedgraph" "3p/\${SAMPLE_ID}_pos3_cpm.bedgraph"
  ln -sf "\${SAMPLE_ID}.3p.neg.cpm.bedgraph" "3p/\${SAMPLE_ID}_neg3_cpm.bedgraph"

  echo "NORMALIZE | LEGACY | Legacy symlinks created"

  ###########################################################################
  # 12) CREATE README
  ###########################################################################

  echo "NORMALIZE | README | Creating documentation..."

  cat > README_normalization.txt <<DOCEOF
================================================================================
NORMALIZED TRACKS — ${sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Normalized coverage tracks for quantitative comparison across samples.
  
  Two normalization methods:
    1. CPM (Counts Per Million): Standard library size normalization
    2. siCPM (Spike-in CPM): Accounts for global changes via spike-in RNA

NORMALIZATION FORMULAS
────────────────────────────────────────────────────────────────────────────

CPM (Counts Per Million):
  CPM = (raw_count / sample_reads) × 1,000,000
  
  Use when:
    • Comparing samples with similar RNA content
    • Standard differential expression analysis
    • No spike-in RNA added

siCPM (Spike-in Normalized CPM):
  siCPM = (raw_count / sample_spike) × (control_spike / control_reads) × 1,000,000
  
  Use when:
    • Global transcription changes expected
    • Spike-in RNA added during library prep
    • Treatment affects overall RNA levels

CONTROL SELECTION (for siCPM)
────────────────────────────────────────────────────────────────────────────
  Condition = "\${CONTROL_LABEL}" (params.control_label), lowest replicate
  number preferred (a merged rep-0 track wins over stray per-replicate rows).

  control_label MUST be set explicitly, per dataset, to your real baseline
  condition name — there is no fallback to "first sample with spike reads"
  if it doesn't match. A non-matching label disables siCPM (factor = 0)
  rather than silently normalizing against the wrong sample.

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample:     ${sample_id}
  Condition:  ${condition}
  Timepoint:  ${timepoint}
  Replicate:  ${replicate}

NORMALIZATION FACTORS
────────────────────────────────────────────────────────────────────────────
  CPM factor:   \${FAC_CPM}
  siCPM factor: \${FAC_SICPM}
  
  \$([ \${SICPM_AVAILABLE} -eq 0 ] && echo "  Note: siCPM disabled (factor = 0)" || echo "  siCPM enabled and available")

FILES
────────────────────────────────────────────────────────────────────────────

Main 3' Tracks (Always Generated):
  3p/${sample_id}.3p.pos.cpm.bedgraph       — Positive strand CPM
  3p/${sample_id}.3p.neg.cpm.bedgraph       — Negative strand CPM
  3p/${sample_id}.3p.pos.cpm.bw             — BigWig format
  3p/${sample_id}.3p.neg.cpm.bw             — BigWig format
  3p/${sample_id}.3p.pos.sicpm.bedgraph     — Positive strand siCPM
  3p/${sample_id}.3p.neg.sicpm.bedgraph     — Negative strand siCPM
  3p/${sample_id}.3p.pos.sicpm.bw           — BigWig format
  3p/${sample_id}.3p.neg.sicpm.bw           — BigWig format

Main 5' Tracks (PE Only):
  5p/${sample_id}.5p.*.cpm.bedgraph         — CPM normalized
  5p/${sample_id}.5p.*.cpm.bw               — BigWig format
  5p/${sample_id}.5p.*.sicpm.bedgraph       — siCPM normalized
  5p/${sample_id}.5p.*.sicpm.bw             — BigWig format

AllMap Tracks (if emit_allmap=true):
  3p/${sample_id}.allMap.3p.*.cpm.bedgraph
  3p/${sample_id}.allMap.3p.*.sicpm.bedgraph
  (Plus corresponding BigWig files)

Legacy Symlinks (for compatibility):
  3p/${sample_id}_pos3_cpm.bedgraph → 3p/${sample_id}.3p.pos.cpm.bedgraph
  3p/${sample_id}_neg3_cpm.bedgraph → 3p/${sample_id}.3p.neg.cpm.bedgraph

Metadata:
  normalization_factors.tsv — CPM and siCPM scaling factors
  tracks_manifest.tsv       — Complete list of all generated tracks
  README_normalization.txt  — This documentation
  normalize_coverage_tracks.log      — Processing log

FILE FORMAT
────────────────────────────────────────────────────────────────────────────
  bedGraph: chr<TAB>start<TAB>end<TAB>normalized_coverage
  BigWig:   Binary indexed format for genome browsers

PROCESSING NOTES
────────────────────────────────────────────────────────────────────────────
  • Single-pass scaling: CPM and siCPM computed together for efficiency
  • No coordinate clipping needed (validated by generate_coverage_tracks module)
  • Negative strand values preserved from upstream mirroring
  • BigWig timeout: \${TIMEOUT_BW} seconds
  • Optional bedGraph sorting: \$([ \${FORCE_SORT} -eq 1 ] && echo "enabled" || echo "disabled")

TRACKS MANIFEST
────────────────────────────────────────────────────────────────────────────
  The tracks_manifest.tsv file contains:
    Column 1: sample       — Sample identifier
    Column 2: end          — 3p or 5p
    Column 3: set          — main or allMap
    Column 4: strand       — pos or neg
    Column 5: scale        — cpm or sicpm
    Column 6: path         — Relative path to bedGraph file

USAGE RECOMMENDATIONS
────────────────────────────────────────────────────────────────────────────

For Genome Browsers:
  • Use BigWig (.bw) files for visualization
  • Load both positive and negative strand tracks
  • Negative values indicate reverse strand coverage

For Downstream Analysis:
  • Use CPM tracks for standard comparisons
  • Use siCPM tracks when spike-in added
  • bedGraph format recommended for computational analysis
  • Manifest file helps programmatic access

For Differential Analysis:
  • CPM tracks suitable for DESeq2, edgeR input
  • siCPM accounts for global expression changes
  • Consider biological replicates for statistics

QUALITY CHECKS
────────────────────────────────────────────────────────────────────────────

Expected Behavior:
  • CPM values typically 0.01 - 1000 for expressed regions
  • siCPM values should be similar to CPM if no global changes
  • All tracks should have same number of intervals as input

Troubleshooting:
  • Empty outputs: Check input bedGraphs exist and non-empty
  • Zero siCPM factor: control_label doesn't match any condition in your
    samplesheet (most common — set it explicitly per dataset), or no
    spike-in reads on this sample/its control. Check the log line
    "NORMALIZE | FACTORS | Cause: ..." above for which one.
  • BigWig conversion failure: Try --force_sort_bedgraph

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  Emit BigWig:      \$([ \${EMIT_BW} -eq 1 ] && echo "Yes" || echo "No")
  Emit siCPM:       \$([ \${EMIT_SICPM} -eq 1 ] && echo "Yes" || echo "No")
  Emit allMap:      \$([ \${EMIT_ALLMAP} -eq 1 ] && echo "Yes" || echo "No")
  Emit 5' tracks:   \$([ \${EMIT_5P} -eq 1 ] && echo "Yes" || echo "No (auto)")
  Control label:    \${CONTROL_LABEL}
  BigWig timeout:   \${TIMEOUT_BW}s
  Force sort:       \$([ \${FORCE_SORT} -eq 1 ] && echo "Yes" || echo "No")

DOWNSTREAM MODULES
────────────────────────────────────────────────────────────────────────────
  These normalized tracks are used by:
  • Divergent transcription detection (allMap 3' tracks)
  • Functional region calling (main 3' tracks)
  • Pol-II metrics calculation (CPM or siCPM depending on spike-in)
  • Quality control and visualization

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: \$(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: ${sample_id}
  Module: 08_normalize_coverage_tracks

================================================================================
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