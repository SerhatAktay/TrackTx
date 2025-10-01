// ============================================================================
// fetch_and_build_index.nf — FASTA fetch + Bowtie2 index (cached)
// ----------------------------------------------------------------------------
// Overview
//   • Deterministic behavior from {genome_id, source, fasta_in, force_rebuild}.
//   • Cross-run reuse via storeDir; results linked via publishDir.
//   • Prefer local ./genomes indices (subdir or flat; e.g., hs1, hs1_primary).
//   • Correct 6-shard checks: 1,2,3,4,rev.1,rev.2 (.bt2 or .bt2l).
//   • Writes REAL files in work dir (no symlink outputs).
//
// Inputs
//   tuple( genome_id, source, fasta_in )
//     • genome_id : UCSC-style assembly (e.g., hg38, mm10, hs1)
//     • source    : note preserved in README (e.g., "UCSC", "custom")
//     • fasta_in  : path to FASTA(.gz) or "-" to force download
//
// Outputs (publishDir):
//   ${params.output_dir}/00_references/${genome_id}/
//     ├── ${genome_id}.fa
//     ├── ${genome_id}.fa.fai
//     ├── ${genome_id}.genome.sizes
//     ├── ${genome_id}.*.bt2*
//     └── README_index.txt
// ============================================================================

nextflow.enable.dsl = 2

process fetch_and_build_index {

  // ── Meta / resources ────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config  
  tag        { genome_id }
  label      'conda'
  cache      'lenient'

  // Ordered, discoverable output location
  storeDir   { "${params.assets_dir ?: "${projectDir}/assets"}/references/${genome_id}" }
  // Publish only lightweight reference artifacts into results; keep heavy .bt2 shards in assets only
  publishDir(
    path: "${params.output_dir}/00_references/${genome_id}",
    mode: 'link',
    saveAs: { filename ->
      def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
      if (name.endsWith('.bt2') || name.contains('.bt2l') || name.matches('.*\\.bt2(\\..*)?$')) return null
      return name
    }
  )

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(genome_id), val(source), path(fasta_in)

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    // Join-friendly metadata tuple
    tuple val(genome_id), val(source), path("${genome_id}.fa"), emit: ref_meta

    // Aux + index shards
    path "${genome_id}.fa.fai"
    path "${genome_id}.genome.sizes"
    path "${genome_id}.*.bt2*", emit: index_files

    // README
    path "README_index.txt"

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  trap 'echo "ERROR [index] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  export LC_ALL=C

  # ── Bindings / locals ────────────────────────────────────────────────────
  ASM="!{genome_id}"
  SRC="!{source}"
  THREADS=!{task.cpus}
  FORCE_REBUILD="!{ params.get('force_rebuild', false) ? 'true' : 'false' }"

  # Container-safe paths: use workflow.workDir or params
  CACHE_ROOT="!{ params.genome_cache ?: '/tmp/genomes_cache' }"
  GENOMES_DIR="${CACHE_ROOT}/local_genomes"
  CACHE_DIR="${CACHE_ROOT}/${ASM}"
  URL_PRIMARY="https://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/bigZips/${ASM}.fa.gz"
  URL_FALLBACK="https://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/bigZips/chromFa.tar.gz"

  mkdir -p "${CACHE_DIR}"
  FASTA_CACHE="${CACHE_DIR}/${ASM}.fa"
  FASTA_DIGEST="${CACHE_DIR}/.${ASM}.fa.sha256"
  IDX_PREFIX="${CACHE_DIR}/${ASM}"

  tmpd=$(mktemp -d -p "${CACHE_DIR}" ".idx_${ASM}.XXXX") || { echo "ERROR  tmpdir" >&2; exit 1; }
  trap 'rm -rf "${tmpd}"' EXIT

  echo "INFO  [index] ▶ asm=${ASM} src=${SRC} cpus=${THREADS} force=${FORCE_REBUILD} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  # ── Tool checks ─────────────────────────────────────────────────────────
  command -v bowtie2-build >/dev/null 2>&1 || { echo "ERROR bowtie2-build not found" >&2; exit 1; }
  command -v samtools     >/dev/null 2>&1 || { echo "ERROR samtools not found"     >&2; exit 1; }
  command -v curl         >/dev/null 2>&1 || { echo "ERROR curl not found"         >&2; exit 1; }
  command -v tar          >/dev/null 2>&1 || { echo "ERROR tar not found"          >&2; exit 1; }

  # ── Helpers ──────────────────────────────────────────────────────────────
  have_full_family() { # $1=prefix, $2=bt2|bt2l
    local p="$1" ext="$2"
    for s in 1 2 3 4 rev.1 rev.2; do [[ -s "${p}.${s}.${ext}" ]] || return 1; done
    return 0
  }

  find_local_prefix(){ # ($basedir, $asm) → best prefix or ""
    local basedir="$1" asm="$2"
    [[ -d "$basedir" ]] || { echo ""; return 0; }
    local cand1="${basedir}/${asm}"
    local cand2="${basedir}/${asm}/${asm}"
    for p in "$cand1" "$cand2"; do
      if have_full_family "$p" bt2 || have_full_family "$p" bt2l; then echo "$p"; return 0; fi
    done
    local f ext pref; shopt -s nullglob
    for ext in bt2 bt2l; do
      for f in "${basedir}/${asm}"*.1."$ext" "${basedir}/${asm}".rev.1."$ext"; do
        [[ -s "$f" ]] || continue
        pref="${f%.1.$ext}"; pref="${pref%.rev}"
        if have_full_family "$pref" "$ext"; then echo "$pref"; shopt -u nullglob; return 0; fi
      done
    done
    shopt -u nullglob
    echo ""
  }

  stage_local_fasta_if_any(){  # NEW: also stage .fai when present
    local pref="$1"
    [[ -s "${FASTA_CACHE}" ]] && return 0
    if   [[ -s "${pref}.fa"    ]]; then
      cp -f "${pref}.fa" "${FASTA_CACHE}"
      [[ -s "${pref}.fa.fai"    ]] && cp -f "${pref}.fa.fai" "${FASTA_CACHE}.fai"   # NEW
    elif [[ -s "${pref}.fa.gz" ]]; then
      gunzip -c "${pref}.fa.gz" > "${FASTA_CACHE}"
    elif [[ -s "${GENOMES_DIR}/${ASM}/${ASM}.fa" ]]; then
      cp -f "${GENOMES_DIR}/${ASM}/${ASM}.fa" "${FASTA_CACHE}"
      [[ -s "${GENOMES_DIR}/${ASM}/${ASM}.fa.fai" ]] && cp -f "${GENOMES_DIR}/${ASM}/${ASM}.fa.fai" "${FASTA_CACHE}.fai"  # NEW
    elif [[ -s "${GENOMES_DIR}/${ASM}/${ASM}.fa.gz" ]]; then
      gunzip -c "${GENOMES_DIR}/${ASM}/${ASM}.fa.gz" > "${FASTA_CACHE}"
    elif [[ -s "${GENOMES_DIR}/${ASM}.fa" ]]; then
      cp -f "${GENOMES_DIR}/${ASM}.fa" "${FASTA_CACHE}"
      [[ -s "${GENOMES_DIR}/${ASM}.fa.fai" ]] && cp -f "${GENOMES_DIR}/${ASM}.fa.fai" "${FASTA_CACHE}.fai"  # NEW
    elif [[ -s "${GENOMES_DIR}/${ASM}.fa.gz" ]]; then
      gunzip -c "${GENOMES_DIR}/${ASM}.fa.gz" > "${FASTA_CACHE}"
    fi
  }

  # ─────────────────────────────────────────────────────────────────────────
  # 1) Fast path — existing cache or local indices
  # ─────────────────────────────────────────────────────────────────────────
  already_built=0
  if [[ "${FORCE_REBUILD}" != "true" ]]; then
    if have_full_family "${IDX_PREFIX}" bt2 || have_full_family "${IDX_PREFIX}" bt2l; then
      echo "INFO  Cache already has full Bowtie2 index → skip build"; already_built=1
    else
      LOCAL_PREFIX="$(find_local_prefix "${GENOMES_DIR}/${ASM}" "${ASM}")"
      [[ -z "${LOCAL_PREFIX}" ]] && LOCAL_PREFIX="$(find_local_prefix "${GENOMES_DIR}" "${ASM}")"
      if [[ -n "${LOCAL_PREFIX}" ]]; then
        echo "INFO  Found local Bowtie2 index at: ${LOCAL_PREFIX} → staging into cache"
        shopt -s nullglob
        for f in "${LOCAL_PREFIX}".*.bt2*; do cp -f "$f" "${CACHE_DIR}/"; done
        shopt -u nullglob
        stage_local_fasta_if_any "${LOCAL_PREFIX}"
        # purge zero-byte shards (partial copies)
        for suf in 1 2 3 4 rev.1 rev.2; do
          for ext in bt2 bt2l; do f="${IDX_PREFIX}.${suf}.${ext}"; [[ -e "$f" && ! -s "$f" ]] && rm -f "$f"; done
        done
        if have_full_family "${IDX_PREFIX}" bt2 || have_full_family "${IDX_PREFIX}" bt2l; then
          echo "INFO  Local indices staged → skip build"; already_built=1
        else
          echo "WARN  Local indices incomplete after staging → will build"
        fi
      else
        echo "INFO  No local indices under ${GENOMES_DIR}"
      fi
    fi
  else
    echo "INFO  force_rebuild=true — ignoring cache/local indices"
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 2) Stage FASTA (input preferred; else UCSC download)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${already_built}" -eq 0 ]]; then
    if [[ -s "${FASTA_CACHE}" && "${FORCE_REBUILD}" != "true" ]]; then
      echo "INFO  Using cached FASTA: ${FASTA_CACHE}"
    else
      if [[ -e "!{fasta_in}" && -s "!{fasta_in}" && "!{fasta_in}" != "-" ]]; then
        echo "INFO  Using provided FASTA input"
        if [[ "!{fasta_in}" == *.gz ]]; then
          gunzip -c -- "!{fasta_in}" > "${tmpd}/${ASM}.fa"
        else
          cp -f -- "!{fasta_in}" "${tmpd}/${ASM}.fa"
        fi
      else
        echo "INFO  Download FASTA from UCSC"
        set +e
        curl -fLsS --retry 5 --retry-delay 3 "${URL_PRIMARY}" | gunzip -c > "${tmpd}/${ASM}.fa"
        rc=$?; set -e
        if [[ $rc -ne 0 || ! -s "${tmpd}/${ASM}.fa" ]]; then
          echo "WARN  Primary failed → fallback chromFa.tar.gz"
          curl -fLsS --retry 5 --retry-delay 3 -o "${tmpd}/chromFa.tar.gz" "${URL_FALLBACK}"
          mkdir -p "${tmpd}/chroms"
          tar -xzf "${tmpd}/chromFa.tar.gz" -C "${tmpd}/chroms"
          cat ${tmpd}/chroms/*.fa > "${tmpd}/${ASM}.fa"
        fi
      fi
      [[ -s "${tmpd}/${ASM}.fa" ]] || { echo "ERROR  FASTA empty" >&2; exit 1; }
      grep -q '^>' "${tmpd}/${ASM}.fa" || { echo "ERROR  Not a FASTA" >&2; exit 1; }
      mv -f "${tmpd}/${ASM}.fa" "${FASTA_CACHE}"; sync
      sha256sum "${FASTA_CACHE}" | awk '{print $1}' > "${FASTA_DIGEST}" || true
    fi
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 2b) Ensure faidx + genome.sizes even on fast path            # NEW BLOCK
  # ─────────────────────────────────────────────────────────────────────────
  if [[ -s "${FASTA_CACHE}" ]]; then
    [[ -s "${FASTA_CACHE}.fai" ]] || samtools faidx "${FASTA_CACHE}"
    cut -f1,2 "${FASTA_CACHE}.fai" > "${CACHE_DIR}/${ASM}.genome.sizes"
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 3) Build Bowtie2 index (size-probe → --large-index; low-mem retry)
  # ─────────────────────────────────────────────────────────────────────────
  if [[ "${already_built}" -eq 0 ]]; then
    size_bytes=$(stat -c %s "${FASTA_CACHE}" 2>/dev/null || stat -f %z "${FASTA_CACHE}")
    [[ "${size_bytes}" -ge $((4*1024*1024*1024)) ]] && LARGE_FLAG="--large-index" || LARGE_FLAG=""

    if ! have_full_family "${IDX_PREFIX}" bt2 && ! have_full_family "${IDX_PREFIX}" bt2l; then
      echo "INFO  Building Bowtie2 index (threads=${THREADS})"
      set +e
      bowtie2-build --threads ${THREADS} ${LARGE_FLAG:+${LARGE_FLAG}} "${FASTA_CACHE}" "${IDX_PREFIX}"
      rc=$?; set -e
      if [[ $rc -ne 0 ]]; then
        echo "WARN  Retry with low-mem flags"
        bowtie2-build --threads 1 ${LARGE_FLAG:+${LARGE_FLAG}} --bmaxdivn 8 --dcv 1024 \
                      "${FASTA_CACHE}" "${IDX_PREFIX}"
      fi
    fi

    have_full_family "${IDX_PREFIX}" bt2 || have_full_family "${IDX_PREFIX}" bt2l \
      || { echo "ERROR  Incomplete index family" >&2; exit 1; }
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 4) Stage declared outputs
  # ─────────────────────────────────────────────────────────────────────────
  # Avoid "cp: not writing through dangling symlink" by unlinking targets first
  rm -f -- "${ASM}.fa" "${ASM}.fa.fai" "${ASM}.genome.sizes" || true
  cp -f -- "${FASTA_CACHE}"                    "${ASM}.fa"
  cp -f -- "${FASTA_CACHE}.fai"                "${ASM}.fa.fai"
  cp -f -- "${CACHE_DIR}/${ASM}.genome.sizes"  "${ASM}.genome.sizes"
  shopt -s nullglob
  for f in "${IDX_PREFIX}".*.bt2*; do cp -f -- "$f" .; done
  shopt -u nullglob

  # ─────────────────────────────────────────────────────────────────────────
  # 5) README
  # ─────────────────────────────────────────────────────────────────────────
  cat > README_index.txt <<TXT
Reference assets for assembly: ${ASM}

Files
  • ${ASM}.fa
  • ${ASM}.fa.fai
  • ${ASM}.genome.sizes
  • ${ASM}.*.bt2*

Notes
  - Local reuse : ./genomes/<ASM>{,/} preferred when present
  - Cache dir   : ${CACHE_DIR}
  - Build flags : ${LARGE_FLAG:+${LARGE_FLAG}} (threads=${THREADS})
  - Force       : ${FORCE_REBUILD}
  - Source      : ${SRC}
TXT

  echo "INFO  [index] ✔ done for ${ASM} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
