#!/usr/bin/env bash
# Run every samplesheet through run_pipeline.sh (conda_server profile).
# Between attempts of one sample, work/ is kept so the retry auto-resumes; it is wiped
# after the sample succeeds or fails MAX_ATTEMPTS times. A failed sample never stops the loop.
set -uo pipefail
cd "$(dirname "$0")" || exit 1
FAILED=()

LOG_FILE=run_all_samples.log
DONE_FILE=run_all_samples.done   # labels = params basename; delete a line to re-run that sample
MAX_ATTEMPTS=2

# "samplesheet params" (basenames, no extension). Comment out (#) to skip.
RUNS="
arabidopsis_leaf_liu2019 arabidopsis_leaf_liu2019
cassava_Nase3_vera2021 cassava_Nase3_vera2021
dog_DH82_himanen2025 dog_DH82_himanen2025
drosophila_S2_duarte2016 drosophila_S2_duarte2016
drosophila_S2Rplus_prajapati2024 drosophila_S2Rplus_prajapati2024
ecoli_MG1655_vill2024 ecoli_MG1655_vill2024
human_K562_dukler2017 human_K562_dukler2017
human_K562_dukler2017 human_K562_dukler2017_hg19
human_K562_vihervaara2017 human_K562_vihervaara2017
human_K562_vihervaara2020 human_K562_vihervaara2020
maize_B73shoot_vera2021 maize_B73shoot_vera2021
mouse_MEF_himanen2022 mouse_MEF_himanen2022
"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

# The results share is CIFS/krb5i: an expired ticket kills open files ("Host is down").
# Renew the TGT in the background. Fails once the ticket's renew-until limit passes.
( while sleep 1800; do kinit -R 2>/dev/null || echo "[$(date '+%F %T')] WARN kinit -R failed" >> "$LOG_FILE"; done ) &
KRB_PID=$!
trap 'kill "$KRB_PID" 2>/dev/null' EXIT
# run_pipeline.sh traps INT itself and exits 130; without this the loop would carry on to the next sample.
trap 'log "Interrupted, stopping."; exit 130' INT TERM

mkdir -p logs_nextflow
[[ -f "$DONE_FILE" ]] || : > "$DONE_FILE"

while read -r sheet name; do
    [[ -z "$sheet" || "$sheet" == \#* ]] && continue
    samplesheet="samplesheets/$sheet.csv"
    params="params/$name.yaml"

    grep -qxF "$name" "$DONE_FILE" && { log "SKIP $name: done"; continue; }
    [[ -f "$samplesheet" ]] || { log "SKIP $name: no $samplesheet"; continue; }
    [[ -f "$params" ]]      || { log "SKIP $name: no $params"; continue; }

    ok=0
    for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
        log "START $name (attempt $attempt/$MAX_ATTEMPTS)"
        ./run_pipeline.sh --samplesheet "$samplesheet" --params-file "$params" \
            -profile conda_server --skip-countdown --no-resume-prompt --no-docker-prompt \
            2>&1 | tee -a "$LOG_FILE"
        status=${PIPESTATUS[0]}
        cp -f .nextflow.log "logs_nextflow/$name.attempt$attempt.log" 2>/dev/null
        [[ $status -eq 130 ]] && { log "Interrupted, stopping."; exit 130; }
        [[ $status -eq 0 ]] && { ok=1; break; }
        log "FAIL $name (exit $status)"
    done

    if [[ $ok -eq 1 ]]; then
        log "OK $name - cleaning work/"
        rm -rf work .nextflow*
        echo "$name" >> "$DONE_FILE"
    else
        # Wipe too: a leftover .nextflow/work would make the NEXT sample auto-resume
        # from this one's state. Logs are kept in logs_nextflow/ and the log file.
        log "GIVING UP on $name after $MAX_ATTEMPTS attempts - continuing with next sample (re-run script to retry it)"
        rm -rf work .nextflow*
        FAILED+=("$name")
    fi
done <<< "$RUNS"

log "Finished; failed: ${FAILED[*]:-none}"
[[ ${#FAILED[@]} -eq 0 ]]
