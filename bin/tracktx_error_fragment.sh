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
