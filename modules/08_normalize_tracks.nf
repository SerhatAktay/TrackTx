// ============================================================================
// 08_normalize_tracks.nf — CPM & siCPM for 3′/5′ (main + allMap), stranded
// ----------------------------------------------------------------------------
// Overview
//   • Computes CPM and (optionally) spike-in CPM (siCPM) per track
//   • Handles main BAM tracks and allMap tracks; 3′ always, 5′ if present
//   • Writes a compact, navigable folder structure + a manifest TSV
//   • Emits legacy 3′ CPM symlinks to avoid breaking downstream wiring
//
// Why it’s faster
//   • Single-pass scaling per input bedGraph → writes CPM and siCPM together
//   • Parallelizes per-track work to utilize available CPUs
//   • Optional sort before BigWig (skip unless forced)
//
// Inputs
//   tuple(
//     sample_id,
//     pos3_bg, neg3_bg,
//     pos5_bg?, neg5_bg?,
//     allmap3p_pos_bg, allmap3p_neg_bg,
//     allmap5p_pos_bg?, allmap5p_neg_bg?,
//     condition, timepoint, replicate,
//     counts_master.tsv,
//     genes   // ignored (iface stability)
//   )
//   genome_fa (FASTA; .fai created if missing)
//
// Outputs (publishDir):
//   ${params.output_dir}/05_normalized_tracks/${sample_id}/
//     ├── 3p/
//     │   ├── main/   /(cpm|sicpm)/(pos|neg).(bedgraph|bw)
//     │   └── allMap/ /(cpm|sicpm)/(pos|neg).(bedgraph|bw)   (if emit_allmap)
//     ├── 5p/ (PE or emit_5p)
//     │   ├── main/   /(cpm|sicpm)/(pos|neg).(bedgraph|bw)
//     │   └── allMap/ /(cpm|sicpm)/(pos|neg).(bedgraph|bw)
//     ├── normalization_factors.tsv
//     ├── tracks_manifest.tsv   # sample,end,set,strand,scale,path
//     ├── README_normalization.txt
//     └── Legacy symlinks (for downstream compatibility):
//         3p/${sample_id}_pos3_cpm.bedgraph → 3p/main/cpm/pos.bedgraph
//         3p/${sample_id}_neg3_cpm.bedgraph → 3p/main/cpm/neg.bedgraph
//
// Feature flags (defaults shown)
//   params.norm.emit_bw      = true
//   params.norm.emit_sicpm   = true
//   params.norm.emit_allmap  = true
//   params.norm.emit_5p      = null  // auto: only if 5′ inputs exist; set true to force
//   params.force_sort_bedgraph = false
//   params.control_label     = 'CTRL'  // for siCPM baseline selection
// ============================================================================

nextflow.enable.dsl = 2

process normalize_tracks {

  // ── Meta / resources ─────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  tag     { sample_id }
  label   'conda'
  cache   'lenient'

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/05_normalized_tracks/${sample_id}",
            mode: 'copy', overwrite: true

  // Tools
  conda (params.conda_norm ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          val(pos3_bg),
          val(neg3_bg),
          val(pos5_bg),
          val(neg5_bg),
          val(am3p_pos_bg),
          val(am3p_neg_bg),
          val(am5p_pos_bg),
          val(am5p_neg_bg),
          val(condition), val(timepoint), val(replicate),
          path(counts_master_tsv),    // Changed from val() to path() to stage file in container
          val(genes_unused)
    path genome_fa

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    // Minimal legacy tuple (kept for downstream wiring)
    tuple val(sample_id),
          path("3p/${sample_id}_pos3_cpm.bedgraph"),
          path("3p/${sample_id}_neg3_cpm.bedgraph"),
          path("normalization_factors.tsv"),
          val(condition), val(timepoint), val(replicate),
          emit: norm_tuple

    // Docs
    path "README_normalization.txt", emit: readme
    path "tracks_manifest.tsv",      emit: manifest

    // Convenience emits (optional; keep light) — updated to flattened paths
    path "3p/${sample_id}.3p.pos.cpm.bedgraph", emit: pos3_cpm_bg
    path "3p/${sample_id}.3p.neg.cpm.bedgraph", emit: neg3_cpm_bg

    // Publish full normalized trees so all expected artifacts are visible
    path "3p/**", emit: tree3p
    path "5p/**", optional: true, emit: tree5p

    // Log (optional)
    path "normalize_tracks.log", optional: true, emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  trap 'echo "ERROR  [normalize] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a normalize_tracks.log) 2>&1

  # ── Bindings / locals ────────────────────────────────────────────────────
  SID='!{sample_id}'
  COND='!{condition}'
  TP='!{timepoint}'
  REP='!{replicate}'
  # counts_master.tsv is now staged as a path input (no path resolution needed)
  COUNTS='!{counts_master_tsv}'
  GENFA='!{genome_fa}'

  POS3='!{pos3_bg}'; NEG3='!{neg3_bg}'
  POS5='!{pos5_bg}'; NEG5='!{neg5_bg}'
  AM3P_POS='!{am3p_pos_bg}'; AM3P_NEG='!{am3p_neg_bg}'
  AM5P_POS='!{am5p_pos_bg}'; AM5P_NEG='!{am5p_neg_bg}'

  # Feature flags (unset or 'true' => enabled)
  if [[ -z "!{ params.norm?.emit_bw }" || "!{ params.norm?.emit_bw }" == "true" ]]; then EMIT_BW=1; else EMIT_BW=0; fi
  if [[ -z "!{ params.norm?.emit_sicpm }" || "!{ params.norm?.emit_sicpm }" == "true" ]]; then EMIT_SICPM=1; else EMIT_SICPM=0; fi
  if [[ -z "!{ params.norm?.emit_allmap }" || "!{ params.norm?.emit_allmap }" == "true" ]]; then EMIT_ALLMAP=1; else EMIT_ALLMAP=0; fi
  FORCE_SORT=$([[ "!{ params.force_sort_bedgraph }" == "true" ]] && echo 1 || echo 0)
  CTRL_LABEL='!{ params.control_label ?: "CTRL" }'
  CLIP_BG=$([[ "!{ params.norm?.clip_bg }" == "true" ]] && echo 1 || echo 0)
  TIMEOUT_BW='!{ params.norm?.timeout_bw ?: 900 }'

  # 5′ emit: auto if unset, else respect user flag
  # In AUTO mode, also emit 5′ if placeholder files exist (size 0) so stubs are created
  EMIT_5P_INIT='!{ params.norm?.emit_5p == null ? "__AUTO__" : (params.norm?.emit_5p in [true, "true"] ? "1" : "0") }'
  if [[ "${EMIT_5P_INIT}" == "__AUTO__" ]]; then
    if [[ -e "${POS5}" || -e "${NEG5}" ]]; then
      EMIT_5P=1
    else
      EMIT_5P=0
    fi
  else
    EMIT_5P=${EMIT_5P_INIT}
  fi

  THREADS=!{task.cpus}

  echo "INFO  [normalize] ▶ sample=${SID} cond=${COND} tp=${TP} rep=${REP} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "INFO  flags: emit_bw=${EMIT_BW} emit_sicpm=${EMIT_SICPM} emit_allmap=${EMIT_ALLMAP} emit_5p=${EMIT_5P} force_sort=${FORCE_SORT}"
  echo "INFO  counts_master: ${COUNTS} (staged from input)"

  [[ -s "${COUNTS}" ]] || { echo "ERROR counts_master.tsv missing or empty: ${COUNTS}" >&2; echo "PWD=$(pwd)" >&2; ls -la "${COUNTS%/*}" >&2 2>/dev/null || true; exit 1; }
  command -v python3 >/dev/null || { echo "ERROR python3 not found" >&2; exit 1; }
  command -v bedGraphToBigWig >/dev/null 2>&1 || [[ "${EMIT_BW}" -eq 0 ]] || { echo "ERROR bedGraphToBigWig not found" >&2; exit 1; }

  mkdir -p 3p 5p

  # ─────────────────────────────────────────────────────────────────────────
  # 1) Compute scaling factors (fac_cpm, fac_sicpm)
  # ─────────────────────────────────────────────────────────────────────────
  python3 - "${COUNTS}" "${SID}" "${CTRL_LABEL}" > factors.tmp <<'PY'
import sys, csv
path, sid, ctrl_label = sys.argv[1], sys.argv[2], sys.argv[3].strip().lower()
rows = list(csv.DictReader(open(path), delimiter='\t'))
def low(x): return (x or '').strip().lower()
def is_ctrl(r):
    return low(r.get('condition')) == ctrl_label and str(r.get('replicate','')).strip() in ('1','r1','R1')
s = next((r for r in rows if r.get('sample')==sid), None)
if not s: print("0\t0"); sys.exit(0)
c_sample  = int(s.get('main_reads',0))
si_sample = int(s.get('spike_reads',0))
ctrl = next((r for r in rows if is_ctrl(r)), None)
if ctrl is None:
    ctrl = next((r for r in rows if int(r.get('spike_reads',0))>0), None)
c_ctrl  = int(ctrl.get('main_reads',0))  if ctrl else 0
si_ctrl = int(ctrl.get('spike_reads',0)) if ctrl else 0
fac_cpm   = (1_000_000.0 / c_sample) if c_sample>0 else 0.0
fac_sicpm = (si_ctrl/float(si_sample))*(1_000_000.0/float(c_ctrl)) if (si_sample>0 and si_ctrl>0 and c_ctrl>0) else 0.0
print(f"{fac_cpm:.10f}\t{fac_sicpm:.10f}")
PY

  read -r FAC_CPM FAC_SICPM < factors.tmp
  [[ "${FAC_CPM}" != "0.0000000000" ]] || { echo "ERROR c_sample==0; cannot compute CPM" >&2; exit 1; }
  if awk -v x="${FAC_SICPM}" 'BEGIN{exit (x>0?0:1)}'; then
    echo "INFO  siCPM enabled: fac_sicpm=${FAC_SICPM}"
  else
    echo "WARN  siCPM disabled (factor=0.0 due to missing/zero baselines)"
  fi
  echo "INFO  factors: fac_cpm=${FAC_CPM} fac_sicpm=${FAC_SICPM}"

  # ─────────────────────────────────────────────────────────────────────────
  # 2) genome.sizes (clip target for BigWig)
  #    Reuse prebuilt index if available next to the source FASTA
  # ─────────────────────────────────────────────────────────────────────────
  GENFA_SRC='!{ genome_fa.toString() }'
  if [[ -s "${GENFA_SRC}.fai" && ! -s "${GENFA}.fai" ]]; then
    ln -sf "${GENFA_SRC}.fai" "${GENFA}.fai" || true
  fi
  if [[ -s "${GENFA}.fai" ]]; then
    cut -f1,2 "${GENFA}.fai" > genome.sizes
  else
    samtools faidx "${GENFA}" >/dev/null 2>&1
    cut -f1,2 "${GENFA}.fai" > genome.sizes
  fi
  LC_ALL=C sort -k1,1 -u -o genome.sizes genome.sizes

  # ─────────────────────────────────────────────────────────────────────────
  # 3) Helpers
  # ─────────────────────────────────────────────────────────────────────────
  clip_bg_inplace() {
    local IN_BG="$1"; local tmp="${1}.cliptmp"
    awk -v OFS='\t' 'FNR==NR{L[$1]=$2; next}
         (NF>=3)&&($0!~/^(track|browser|#)/){
           s=$2+0; e=$3+0; if(s<0)s=0; m=L[$1]+0; if(m==0) next;
           if(e>m)e=m; if(s<e){ if(NF>=4) print $1,s,e,$4; else print $1,s,e; }
         }' genome.sizes "$IN_BG" > "$tmp" && mv -f "$tmp" "$IN_BG"
  }

  make_bw () {
    local IN_BG="$1" BW="$2"
    [[ -s "$IN_BG" ]] || { : > "$BW"; return 0; }
    echo "INFO  [normalize] make_bw: preparing $(basename \"$IN_BG\") → $(basename \"$BW\")"
    if [[ ${EMIT_BW} -eq 1 && ${CLIP_BG} -eq 1 ]]; then
      echo "INFO  [normalize] make_bw: clipping to genome bounds"
      clip_bg_inplace "$IN_BG"
    else
      echo "INFO  [normalize] make_bw: skipping clip (emit_bw=${EMIT_BW} clip_bg=${CLIP_BG})"
    fi
    if [[ ${FORCE_SORT} -eq 1 ]]; then
      echo "INFO  [normalize] make_bw: sorting bedGraph prior to BigWig"
      LC_ALL=C sort -k1,1 -k2,2n -o "$IN_BG" "$IN_BG"
    fi
    if [[ ${EMIT_BW} -eq 1 ]]; then
      local LINES=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "$IN_BG" 2>/dev/null || echo 0)
      echo "INFO  [normalize] make_bw: bedGraphToBigWig (rows=${LINES}, timeout=${TIMEOUT_BW}s)"
      if timeout "${TIMEOUT_BW}" bedGraphToBigWig "$IN_BG" genome.sizes "$BW"; then
        local BWSZ=$(stat -f%z "$BW" 2>/dev/null || stat -c%s "$BW" 2>/dev/null || echo "unknown")
        echo "INFO  [normalize] make_bw: wrote $(basename \"$BW\") (${BWSZ} bytes)"
      else
        echo "WARN  [normalize] make_bw: bedGraphToBigWig failed/timed out; emitting empty $(basename \"$BW\")"
        : > "$BW"
      fi
    else
      : > "$BW"
    fi
  }

  # Single-pass scaler: reads BG once, writes CPM and (optionally) siCPM
  scale_bg () {
    # $1=in_bg  $2=end(3p|5p)  $3=pos|neg  $4=set(main|allMap)
    local IN="$1" END="$2" STR="$3" SETLBL="$4"
    local SETSUF=""
    if [[ "$SETLBL" == "allMap" ]]; then SETSUF=".allMap"; fi
    local PREFIX="${SID}${SETSUF}.${END}.${STR}"
    local OUT_CPM_BG="${END}/${PREFIX}.cpm.bedgraph"
    local OUT_CPM_BW="${END}/${PREFIX}.cpm.bw"
    local OUT_SI_BG="${END}/${PREFIX}.sicpm.bedgraph"
    local OUT_SI_BW="${END}/${PREFIX}.sicpm.bw"

    if [[ -s "$IN" ]]; then
      echo "INFO  [normalize] scale: ${END} ${SETLBL} ${STR} — input=$(basename \"$IN\")"
      local IN_ROWS=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "$IN" 2>/dev/null || echo 0)
      echo "INFO  [normalize] scale: input rows=${IN_ROWS}"
      local T0=$(date +%s)
      if [[ ${EMIT_SICPM} -eq 1 && $(awk -v x="${FAC_SICPM}" 'BEGIN{print (x>0)?1:0}') -eq 1 ]]; then
        echo "INFO  [normalize] scale: writing CPM and siCPM (fac_cpm=${FAC_CPM}, fac_sicpm=${FAC_SICPM})"
        awk -v fc="${FAC_CPM}" -v fs="${FAC_SICPM}" -v OFS='\t' '
          (NF>=4)&&($0!~/^(track|browser|#)/){
            c=$4*fc; s=$4*fs;
            print $1,$2,$3,c > "'""$OUT_CPM_BG""'";
            print $1,$2,$3,s > "'""$OUT_SI_BG""'";
          }' "$IN"
        echo "INFO  [normalize] scale: wrote $(basename \"$OUT_CPM_BG\"), $(basename \"$OUT_SI_BG\")"
      else
        echo "INFO  [normalize] scale: writing CPM only (fac_cpm=${FAC_CPM}); siCPM disabled"
        awk -v fc="${FAC_CPM}" -v OFS='\t' '
          (NF>=4)&&($0!~/^(track|browser|#)/){ print $1,$2,$3,$4*fc }' "$IN" > "$OUT_CPM_BG"
        : > "$OUT_SI_BG"
        echo "INFO  [normalize] scale: wrote $(basename \"$OUT_CPM_BG\")"
      fi
      local T1=$(date +%s)
      local CPM_ROWS=$(awk 'BEGIN{n=0} $0!~/^(track|browser|#)/{n++} END{print n}' "$OUT_CPM_BG" 2>/dev/null || echo 0)
      local DT=$(( T1 - T0 ))
      if [[ $DT -lt 1 ]]; then DT=1; fi
      local RPS=$(( IN_ROWS / DT ))
      echo "INFO  [normalize] scale: normalize done in ${DT}s (~${RPS} rows/s)"

      make_bw "$OUT_CPM_BG" "$OUT_CPM_BW"
      [[ -e "$OUT_SI_BG" ]] && make_bw "$OUT_SI_BG" "$OUT_SI_BW"

      # Manifest rows (flattened paths)
      echo -e "${SID}\t${END}\t${SETLBL}\t${STR}\tcpm\t${OUT_CPM_BG}" >> tracks_manifest.tsv
      [[ -s "$OUT_SI_BG" ]] && echo -e "${SID}\t${END}\t${SETLBL}\t${STR}\tsicpm\t${OUT_SI_BG}" >> tracks_manifest.tsv
    else
      : > "$OUT_CPM_BG"; : > "$OUT_CPM_BW"; : > "$OUT_SI_BG"; : > "$OUT_SI_BW"
    fi
    echo "INFO  [normalize] scale: ${END} ${SETLBL} ${STR} — done"
  }

  # ─────────────────────────────────────────────────────────────────────────
  # 4) Normalize main 3′/5′ and allMap 3′/5′ (parallel where safe)
  # ─────────────────────────────────────────────────────────────────────────
  : > tracks_manifest.tsv

  echo "INFO  [normalize] 3′ main"
  scale_bg "${POS3}" "3p"  "pos" "main"
  scale_bg "${NEG3}" "3p"  "neg" "main"
  echo "INFO  [normalize] 3′ main done"

  if [[ ${EMIT_5P} -eq 1 ]]; then
    echo "INFO  [normalize] 5′ main"
    scale_bg "${POS5}" "5p"  "pos" "main"
    scale_bg "${NEG5}" "5p"  "neg" "main"
    echo "INFO  [normalize] 5′ main done"
  fi

  if [[ ${EMIT_ALLMAP} -eq 1 ]]; then
    echo "INFO  [normalize] 3′ allMap"
    scale_bg "${AM3P_POS}" "3p" "pos" "allMap"
    scale_bg "${AM3P_NEG}" "3p" "neg" "allMap"
    echo "INFO  [normalize] 3′ allMap done"

    if [[ ${EMIT_5P} -eq 1 ]]; then
      echo "INFO  [normalize] 5′ allMap"
      scale_bg "${AM5P_POS}" "5p" "pos" "allMap"
      scale_bg "${AM5P_NEG}" "5p" "neg" "allMap"
      echo "INFO  [normalize] 5′ allMap done"
    fi
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 5) Factors, README, and legacy symlinks for 3′ CPM (pos/neg)
  # ─────────────────────────────────────────────────────────────────────────
  printf "method\tfactor\nCPM\t%s\nsiCPM\t%s\n" "${FAC_CPM}" "${FAC_SICPM}" > normalization_factors.tsv

  # Legacy symlinks for compatibility with downstream emits
  mkdir -p 3p
  ln -sf "${SID}.3p.pos.cpm.bedgraph" "3p/${SID}_pos3_cpm.bedgraph"
  ln -sf "${SID}.3p.neg.cpm.bedgraph" "3p/${SID}_neg3_cpm.bedgraph"

  cat > README_normalization.txt <<TXT
Normalization for sample: ${SID}

Formulas
  CPM   = score * (1e6 / c_sample)
  siCPM = score * (si_ctrl / si_sample) * (1e6 / c_ctrl)

Control selection for siCPM
  Preferred: condition == ${CTRL_LABEL} & replicate == 1
  Fallback : first row in counts_master.tsv with spike_reads > 0

Notes
  • NEG tracks are already mirrored upstream; scaling preserves sign.
  • 3′ is always produced; 5′ only if inputs present (or emit_5p=true).
  • BigWigs built via bedGraphToBigWig after bounds clipping; set
    params.force_sort_bedgraph=true if your system complains about sort order.
  • Manifest "tracks_manifest.tsv" lists every produced bedGraph (and scale).
TXT

  echo "INFO  [normalize] ✔ ${SID} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
