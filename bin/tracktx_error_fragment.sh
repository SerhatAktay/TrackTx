#!/usr/bin/env bash
# =============================================================================
# TrackTx Error Reporting — Template for pipeline modules
# =============================================================================
# Copy this function into process scripts for consistent, readable error output.
# Errors go to stderr so they appear clearly in Nextflow's "Command error" section.
#
# Usage: tracktx_error "module_name" "problem description" "fix instruction" [exit_code]
#
# Example:
#   tracktx_error "detect_divergent_transcription" "Missing Python deps" "pip install numpy pandas"
#
# Use log-only stdout to reduce Nextflow error wall of text on failure:
#   exec > module.log
#   exec 2> >(tee -a module.log >&2)
# Add ERR trap for unexpected failures:
#   trap 'tracktx_error "module_name" "Unexpected process failure" "Check module.log in work dir"' ERR
# =============================================================================

tracktx_error() {
  local module="$1" problem="$2" fix="$3" code="${4:-1}"
  echo "" >&2
  echo "═══════════════════════════════════════════════════════════════════════" >&2
  echo "TRACKTX ERROR" >&2
  echo "═══════════════════════════════════════════════════════════════════════" >&2
  echo "Module:  ${module}" >&2
  echo "Problem: ${problem}" >&2
  echo "Fix:     ${fix}" >&2
  echo "═══════════════════════════════════════════════════════════════════════" >&2
  exit "$code"
}

# =============================================================================
# tracktx_size() — cross-platform file size (GNU stat vs BSD/macOS stat)
# =============================================================================
# Replaces the `stat -c%s "$f" 2>/dev/null || stat -f%z "$f" 2>/dev/null ||
# echo "unknown"` one-liner previously copy-pasted ~57 times across modules.
# Usage: SIZE=$(tracktx_size "${FILE}")            # fallback: "unknown"
#        SIZE=$(tracktx_size "${FILE}" 0)           # fallback: 0 (numeric context)
tracktx_size() {
  stat -c%s "$1" 2>/dev/null || stat -f%z "$1" 2>/dev/null || echo "${2:-unknown}"
}

# =============================================================================
# tracktx_resolve_python() — find a Python 3 interpreter across execution
# environments (Docker/Singularity container, conda profile, bare PATH)
# =============================================================================
# Replaces the identical micromamba/opt-conda/python3 fallback chain
# previously copy-pasted across 8 modules. Sets the global PYTHON_CMD.
# Usage: tracktx_resolve_python
tracktx_resolve_python() {
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  else
    PYTHON_CMD="python3"
  fi
}

# =============================================================================
# tracktx_io_lock_init() / with_io_lock() — serialize whole-BAM I/O passes
# across concurrently running tasks sharing a slow (USB/HDD/NFS) work volume.
# =============================================================================
# Replaces the identical IO_LOCK_FILE/with_io_lock() block previously
# copy-pasted across modules 05, 05b, 06, 07, 13. See generate_coverage_tracks.nf
# for the full writeup of why this lock exists. Usage:
#   tracktx_io_lock_init "${projectDir}/.tracktx_io.lock"
#   with_io_lock samtools index -@ "${THREADS}" "${BAM}"
tracktx_io_lock_init() {
  IO_LOCK_FILE="$1"
  IO_LOCK_TIMEOUT="${TRACKS_IO_LOCK_TIMEOUT:-1800}"
}

# Lock via fd in a subshell, not `flock FILE CMD`: the latter execs CMD and
# cannot run shell functions (e.g. tracktx_stage_immutable). A lock timeout
# exits the subshell non-zero, so the caller's ERR trap still fires.
with_io_lock() {
  (
    flock -w "${IO_LOCK_TIMEOUT}" 9 || exit 1
    "$@"
  ) 9>>"${IO_LOCK_FILE}"
}

# =============================================================================
# tracktx_stage_immutable() — materialize a large, never-mutated-in-place
# source file at a new path as cheaply as the filesystem allows.
# =============================================================================
# Hardlinks (instant, no extra disk) when source and destination share a
# device; falls back to a real copy when that fails for any reason (cross-
# device staging, or a filesystem without hardlink support, e.g. exFAT, or a
# multi-volume network mount where EXDEV is common). Never symlinks: some
# destinations (SMB/CIFS publish targets) reject the target being a symlink
# outright. Only use this for files the caller will never edit in place --
# a hardlink shares the same inode as the source.
# Usage: tracktx_stage_immutable "$SRC" "$DST"
tracktx_stage_immutable() {
  local src="$1" dst="$2"
  ln -f "${src}" "${dst}" 2>/dev/null || cp -f "${src}" "${dst}"
}
