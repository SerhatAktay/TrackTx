// ============================================================================
// capture_tool_versions.nf — Pipeline Provenance / Tool Versions Manifest
// ============================================================================
//
// Purpose:
//   Records the exact tool/package versions used for a run, so a paper's
//   Methods section (or anyone rerunning the pipeline later) has a single
//   authoritative file to cite instead of grepping .nextflow.log.
//
// Why one snapshot, not per-process versions.yml (the nf-core convention):
//   TrackTx deliberately uses a single container/conda environment
//   (envs/Dockerfile, envs/tracktx.yaml) for every process -- see README
//   "Why one container for every step?". Every task in a run therefore
//   shares the exact same toolchain, so one environment dump here is
//   equivalent to capturing it per-process, without the redundancy.
//
// Inputs:
//   None (driven by params + the environment this task runs in)
//
// Outputs:
//   ${params.output_dir}/pipeline_info/versions.yml
//
// ============================================================================

process capture_tool_versions {

  tag        'versions'
  label      'lightweight'
  cache      'lenient'
  conda      "${projectDir}/envs/tracktx.yaml"

  publishDir "${params.output_dir}/pipeline_info",
             mode: params.publish_mode,
             overwrite: true

  output:
    path("versions.yml"), emit: versions

  script:
  """
  #!/usr/bin/env bash
  set -uo pipefail

  {
    echo "# TrackTx pipeline provenance -- generated $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
    echo "pipeline:"
    echo "  name: tracktx"
    echo "  version: '${workflow.manifest.version}'"
    echo "  nextflow_version: '${nextflow.version}'"
    echo "  run_name: '${workflow.runName}'"
    echo "  command_line: '${workflow.commandLine.replace("'", "\\'")}'"
    echo "environment:"

    # Report exact resolved package versions for the ONE environment every
    # process in this run shares (see header comment above). Falls through
    # micromamba (container image) -> conda (conda profile) -> pip, so this
    # works regardless of execution profile.
    if command -v micromamba >/dev/null 2>&1; then
      echo "  manager: micromamba"
      echo "  packages: |"
      micromamba list 2>/dev/null | sed 's/^/    /'
    elif command -v conda >/dev/null 2>&1; then
      echo "  manager: conda"
      echo "  packages: |"
      conda list 2>/dev/null | sed 's/^/    /'
    else
      echo "  manager: pip"
      echo "  packages: |"
      pip list 2>/dev/null | sed 's/^/    /'
    fi
  } > versions.yml

  echo "VERSIONS | wrote \$(wc -l < versions.yml) lines to versions.yml"
  """
}
