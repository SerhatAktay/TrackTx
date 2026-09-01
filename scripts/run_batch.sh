#!/usr/bin/env bash
#
# run_batch.sh -- Sequential multi-sample TrackTx batch runner
#
# PURPOSE
#   Runs ./run_pipeline.sh once per (samplesheet, params-file) pair listed
#   in a jobs file, one sample after another. After each run it:
#     1. deletes the Nextflow work directory and Nextflow's own state
#        (work/, .nextflow/, .nextflow.log*, .nxf_temp/) so every sample
#        starts clean and disk space is reclaimed between (large) PRO-seq
#        runs;
#     2. restarts Docker Desktop and BLOCKS until the daemon is confirmed
#        ready -- not just "the process exists", but that `docker info`
#        responds AND a throwaway container actually runs -- before
#        starting the next sample.
#
#   The container-run check in step 2 is the fix for the previous version
#   of this workflow: it restarted Docker but moved on before Docker had
#   actually finished coming back up (the CLI can respond to `docker info`
#   a few seconds before the VM's file-sharing is usable), so the next
#   sample's first task failed and the batch silently stalled.
#
# USAGE
#   ./scripts/run_batch.sh [jobs_file]
#   (jobs_file defaults to scripts/batch_jobs.txt)
#
#   Run from the pipeline root (the directory containing run_pipeline.sh):
#     cd /path/to/tracktx && ./scripts/run_batch.sh
#
#   For a multi-hour/multi-day batch, run it under nohup/tmux/screen so it
#   survives a closed terminal:
#     nohup ./scripts/run_batch.sh > /dev/null 2>&1 &
#   (per-sample output already goes to logs/run_batch/<timestamp>/ -- see
#   OUTPUT below -- so stdout here is just the same lines duplicated.)
#
# INPUT: jobs file
#   Plain text, one job per line:  <samplesheet path>  <params-file path>
#   Paths are relative to the pipeline root, exactly as you'd pass them to
#   run_pipeline.sh directly. Blank lines and lines starting with # are
#   ignored. See scripts/batch_jobs.example.txt for a template listing the
#   samplesheet/params pairs available at the time it was written.
#
# BEHAVIOUR ON FAILURE
#   If a sample's pipeline run fails, the failure is logged to the summary
#   and the batch moves on to the next sample anyway (cleaning up state and
#   restarting Docker first, same as on success). This means a failed run's
#   work/ and Nextflow log are deleted along with everything else -- the
#   per-sample log file captured below is your only record, so check it
#   before re-running that sample.
#
# OUTPUT
#   - logs/run_batch/<timestamp>/<n>_<sample>.log  (full run_pipeline.sh
#     output for sample n)
#   - logs/run_batch/<timestamp>/summary.tsv (one row per sample: status,
#     exit code, start/end time, log path)
#   - A short summary printed at the end
#
# REQUIREMENTS
#   - macOS with Docker Desktop installed at /Applications/Docker.app.
#     restart_docker() below is macOS-specific (uses osascript/open -a) --
#     see that function to adapt for Linux, where Docker usually runs as a
#     system service (`systemctl restart docker`) instead.
#   - Run from the pipeline root; run_pipeline.sh must be present there.
#   - No other Nextflow run may already be active in this directory -- the
#     script refuses to start if `pgrep` finds one (see the pre-flight
#     check below), since sharing one work/ dir between two live runs
#     corrupts both.
#
# NOTES
#   - EXTRA_PIPELINE_ARGS defaults to (--resume) below: every sample is run
#     with -resume. For a sample with no prior partial run this is a
#     harmless no-op (Nextflow just starts fresh); for a sample that was
#     interrupted mid-run it continues from the last completed task instead
#     of redoing it. Remove it from EXTRA_PIPELINE_ARGS once you no longer
#     need it, or if you'd rather have jobs that redo everything.
#
# No credentials, machine-specific paths, or personal data are used here --
# every path is relative to the pipeline root you run this from.

set -uo pipefail   # deliberately not -e: a failing pipeline run must be
                    # captured and handled, not abort this script.

# ── Configuration (edit as needed) ──────────────────────────────────────
JOBS_FILE="${1:-scripts/batch_jobs.txt}"
DOCKER_APP_NAME="Docker"          # `open -a "$DOCKER_APP_NAME"`
DOCKER_READY_TIMEOUT=600          # seconds to wait for Docker after a restart
DOCKER_POLL_INTERVAL=5            # seconds between readiness checks
DOCKER_QUIT_TIMEOUT=60            # seconds to wait for a clean quit before forcing
EXTRA_PIPELINE_ARGS=(--resume)    # applied to every job -- see NOTES above.
                                   # Add more flags the same way, e.g.
                                   # EXTRA_PIPELINE_ARGS=(--resume --external-drive)

PIPELINE_ROOT="$(pwd)"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG_DIR="logs/run_batch/${TIMESTAMP}"
SUMMARY_FILE="${RUN_LOG_DIR}/summary.tsv"

# ── Small helpers ────────────────────────────────────────────────────────
ts()   { date '+%Y-%m-%d %H:%M:%S'; }
log()  { echo "[$(ts)] $*"; }
warn() { echo "[$(ts)] WARNING: $*" >&2; }
err()  { echo "[$(ts)] ERROR: $*" >&2; }

# ── Sanity checks ────────────────────────────────────────────────────────
if [[ ! -x "./run_pipeline.sh" ]]; then
    err "run_pipeline.sh not found (or not executable) in $(pwd)."
    err "Run this script from the TrackTx pipeline root, e.g.:"
    err "  cd /path/to/tracktx && ./scripts/run_batch.sh"
    exit 1
fi

if [[ ! -f "$JOBS_FILE" ]]; then
    err "Jobs file not found: $JOBS_FILE"
    err "Copy scripts/batch_jobs.example.txt to $JOBS_FILE and edit it, or pass a path:"
    err "  ./scripts/run_batch.sh path/to/jobs.txt"
    exit 1
fi

TOTAL_JOBS=$(grep -Ecv '^[[:space:]]*(#|$)' "$JOBS_FILE")
if [[ "$TOTAL_JOBS" -eq 0 ]]; then
    err "No jobs found in $JOBS_FILE (only blank/comment lines)."
    exit 1
fi

# Refuse to start if a Nextflow run is already active anywhere on this
# machine. Two Nextflow sessions sharing one work/ dir will corrupt both --
# this guard exists because that's almost exactly what nearly happened
# during development of this script.
if pgrep -f "nextflow run" >/dev/null 2>&1; then
    err "A Nextflow process is already running (found via 'pgrep -f \"nextflow run\"')."
    err "Refusing to start: this batch would delete work/ and .nextflow/ out from"
    err "under it the moment its own sample finishes. Let it finish (or stop it"
    err "deliberately) before running this script."
    exit 1
fi

mkdir -p "$RUN_LOG_DIR"
printf 'n\tsamplesheet\tparams_file\tstatus\texit_code\tstarted\tfinished\tlog\n' > "$SUMMARY_FILE"

# ── Docker helpers (macOS / Docker Desktop) ─────────────────────────────
docker_ready() {
    # `docker info` responding is not enough on its own -- it can succeed
    # a few seconds before the VM's file-sharing is actually usable, which
    # is exactly what caused the previous script's batch to stall. Actually
    # running a throwaway container confirms Docker is truly ready.
    docker info >/dev/null 2>&1 || return 1
    docker run --rm hello-world >/dev/null 2>&1
}

wait_for_docker() {
    local waited=0
    log "Waiting for Docker to be ready (timeout ${DOCKER_READY_TIMEOUT}s)..."
    while (( waited < DOCKER_READY_TIMEOUT )); do
        if docker_ready; then
            log "Docker is ready after ${waited}s."
            return 0
        fi
        sleep "$DOCKER_POLL_INTERVAL"
        waited=$(( waited + DOCKER_POLL_INTERVAL ))
    done
    err "Docker did not become ready within ${DOCKER_READY_TIMEOUT}s."
    return 1
}

restart_docker() {
    log "Restarting Docker Desktop..."
    osascript -e "quit app \"${DOCKER_APP_NAME}\"" >/dev/null 2>&1 || true

    local waited=0
    while pgrep -f "Docker.app" >/dev/null 2>&1; do
        sleep 2
        waited=$(( waited + 2 ))
        if (( waited >= DOCKER_QUIT_TIMEOUT )); then
            warn "Docker didn't quit within ${DOCKER_QUIT_TIMEOUT}s, forcing quit..."
            pkill -f "Docker.app" >/dev/null 2>&1 || true
            break
        fi
    done

    sleep 3
    open -a "$DOCKER_APP_NAME"
    wait_for_docker
}

# ── Cleanup between samples ─────────────────────────────────────────────
clean_pipeline_state() {
    cd "$PIPELINE_ROOT" || { err "Lost pipeline root directory"; exit 1; }
    log "Cleaning work/ and Nextflow state..."
    rm -rf work .nextflow .nxf_temp
    rm -f .nextflow.log .nextflow.log.*
}

# ── Pre-flight: make sure Docker is up before the first sample ─────────
log "Checking Docker before starting the batch..."
if docker_ready; then
    log "Docker already ready."
elif ! restart_docker; then
    err "Docker is not available and could not be started. Aborting."
    exit 1
fi

# ── Main loop ────────────────────────────────────────────────────────────
n=0
n_ok=0
n_failed=0

while IFS=$' \t' read -r samplesheet params_file _rest; do
    [[ -z "$samplesheet" || "$samplesheet" == \#* ]] && continue
    n=$(( n + 1 ))

    base="$(basename "$samplesheet" .csv)_$(basename "$params_file" .yaml)"
    sample_log="${RUN_LOG_DIR}/${n}_${base}.log"
    started="$(ts)"

    log "──────────────────────────────────────────────────────────────"
    log "[$n/$TOTAL_JOBS] Starting: $samplesheet + $params_file"
    log "    Log: $sample_log"

    cd "$PIPELINE_ROOT" || exit 1
    ./run_pipeline.sh \
        --samplesheet "$samplesheet" \
        --params-file "$params_file" \
        --no-clear --resume --no-docker-prompt \
        "${EXTRA_PIPELINE_ARGS[@]}" \
        < /dev/null > "$sample_log" 2>&1
    exit_code=$?
    finished="$(ts)"

    if [[ $exit_code -eq 0 ]]; then
        status="OK"
        n_ok=$(( n_ok + 1 ))
        log "[$n/$TOTAL_JOBS] Completed OK: $samplesheet"
    else
        status="FAILED"
        n_failed=$(( n_failed + 1 ))
        err "[$n/$TOTAL_JOBS] Pipeline FAILED (exit $exit_code): $samplesheet -- see $sample_log"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$n" "$samplesheet" "$params_file" "$status" "$exit_code" "$started" "$finished" "$sample_log" \
        >> "$SUMMARY_FILE"

    clean_pipeline_state

    if (( n < TOTAL_JOBS )); then
        if ! restart_docker; then
            err "Docker never came back up after sample $n. Stopping the batch here"
            err "(remaining jobs were NOT attempted) -- fix Docker, then rerun with the"
            err "remaining lines of $JOBS_FILE."
            break
        fi
    else
        log "Last job in the batch -- skipping the final Docker restart."
    fi

done < "$JOBS_FILE"

log "──────────────────────────────────────────────────────────────"
log "Batch finished: $n_ok OK, $n_failed FAILED, out of $n attempted (of $TOTAL_JOBS listed)."
log "Summary: $SUMMARY_FILE"
[[ $n_failed -gt 0 ]] && exit 1
exit 0
