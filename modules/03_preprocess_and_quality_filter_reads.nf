// ============================================================================
// preprocess_and_quality_filter_reads.nf — FASTQ Preprocessing and Quality Control
// ============================================================================
//
// Purpose:
//   Single-pass preprocessing of raw FASTQ files for PRO-seq analysis
//
// Features:
//   • One-pass cutadapt: adapters, barcodes, and length filtering
//   • Barcode/UMI extraction on EITHER read (R1 or R2) and EITHER end (5'/3')
//   • Optional QC-based detection of barcode/UMI read+end+length, with two
//     modes layered on top of the default of trusting your own settings:
//       - "auto"   : "I don't know" — ignore length/location, let QC decide
//       - "verify" : "I'm not sure" — use my settings, but HALT if the QC
//                    data clearly disagrees with them
//   • Quality control on raw and final reads (FastQC)
//   • Comprehensive statistics and reporting
//   • Fast execution with multi-threading
//
// Why Single-Pass?
//   Traditional multi-pass approach:
//     raw → adapter trim → barcode trim → UMI extract → final
//   
//   Our optimized approach:
//     raw → cutadapt (adapters + barcodes + minlen) → UMI extract → final
//   
//   Benefits: Faster execution, fewer intermediate files, simpler logic
//   (detect_mode "trust", the default, keeps this true — detection is
//   skipped entirely unless a slot asks for "auto" or "verify".)
//
// Workflow:
//   1. Optional QC on raw reads (FastQC)
//   2. Optional barcode/UMI QC detection (only if any detect_mode != "trust")
//   3. Single-pass cutadapt (adapters + barcodes + length filter)
//   4. UMI extraction (if enabled)
//   5. QC on final cleaned reads (FastQC)
//   6. Generate statistics and reports
//
// Inputs:
//   tuple(sample_id, reads, condition, timepoint, replicate)
//   val(data_type) : "SE" or "PE"
//
// Outputs:
//   ${params.output_dir}/01_trimmed_fastq/${sample_id}/
//     ├── final_R1.fastq                 — Cleaned R1
//     ├── final_R2.fastq                 — Cleaned R2 (or stub for SE)
//     ├── fastqc_raw/                    — Raw QC reports
//     ├── fastqc_final/                  — Final QC reports
//     ├── trim_stats.tsv                 — Trimming statistics
//     ├── umi_stats.tsv                  — UMI extraction stats (if enabled)
//     ├── barcode_umi_profile.json       — QC detection profile (if any detect_mode != trust)
//     ├── barcode_umi_detect.log         — Detection decisions log (same condition)
//     ├── cutadapt_report.txt            — Detailed cutadapt log
//     ├── umi_extract.log                — UMI extraction log (if enabled)
//     ├── preprocess_reads.log           — Complete process log
//     └── README_01_trimmed_fastq.txt    — Documentation
//
// Parameters:
//   params.fastqc_raw                     : Run QC on raw reads (default: true)
//   params.qc.enabled                     : Enable QC (default: true)
//   params.adapter_trimming.enabled       : Enable adapter trimming
//   params.adapter_trimming.adapter1/2    : Adapter sequences
//   params.adapter_trimming.minlen        : Minimum length (user override)
//
//   params.barcode.enabled                : Enable primary barcode removal
//   params.barcode.read                   : Which read: "R1" or "R2" (default: "R1")
//   params.barcode.location               : Which end: "5" or "3" (default: "5")
//   params.barcode.length                 : Barcode length in bp
//   params.barcode.detect_mode            : "trust" (default) | "auto" | "verify"
//   params.barcode.detect_tolerance       : +/- bp allowed before "verify" halts (default: 1)
//   params.barcode.enabled2/read2/location2/length2 : optional SECOND barcode
//     (e.g. a different barcode on R1 vs R2 in the same library) — same
//     read/location/length/detect_mode model as the primary one, just with
//     a "2" suffix. detect_mode2 controls it independently.
//   params.barcode.detect_mode2 / detect_tolerance2 : as above, for barcode2
//
//   params.umi.enabled                    : Enable UMI extraction
//   params.umi.read                       : Which read: "R1" or "R2" (default: "R1")
//   params.umi.location                   : Which end: "5" or "3" (default: "5")
//   params.umi.length                     : UMI length in bp
//   params.umi.detect_mode                : "trust" (default) | "auto" | "verify"
//   params.umi.detect_tolerance           : +/- bp allowed before "verify" halts (default: 1)
//
//   params.advanced.prep_insert_minlen    : Final insert minimum (default: 12)
//
// detect_mode, in plain terms (set independently per barcode / barcode2 / umi):
//   "trust"  : use read/location/length exactly as given. No QC detection
//              runs for that slot. This is the default, and matches every
//              dataset built before this feature existed — nothing changes
//              for them unless a params.yaml opts in.
//   "auto"   : "I don't know where it is." Ignore read/location/length and
//              let bin/detect_barcode_umi.py's scan of the raw reads propose
//              them instead. If it finds nothing conclusive, that slot is
//              disabled (with a clear warning in the log) rather than
//              guessing.
//   "verify" : "I think I know, but I'm not fully sure." Use read/location/
//              length exactly as given (identical trimming to "trust"), but
//              also run the QC scan and compare it to the claim. If they
//              disagree by more than detect_tolerance bp, the process HALTS
//              via tracktx_error instead of silently proceeding on a claim
//              the data doesn't support. On a match, behaves like "trust".
//
// ============================================================================


process preprocess_and_quality_filter_reads {

  tag    { sample_id }
  label  'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/01_trimmed_fastq/${sample_id}" },
             mode: params.publish_mode,
             overwrite: true,
             saveAs: { filename ->
               if (params.get('publish_trimmed_fastq')?.toString() == 'false') return null  // Skip entire folder to save ~900 MB/sample
               def name = filename instanceof Path ? filename.getFileName().toString() : filename.toString()
               // Only publish processed outputs, NOT the raw input FASTQ files
               // This prevents duplicating large raw FASTQ files in the results folder
               // Publish only: final_*.fastq, QC reports, stats, logs, and README
               if (name.startsWith('final_') ||
                   name.startsWith('fastqc_') ||
                   name.contains('_stats.tsv') ||
                   name.contains('_report.txt') ||
                   name.contains('_extract.log') ||
                   name.contains('preprocess_reads.log') ||
                   name.startsWith('barcode_umi_') ||
                   name.startsWith('README_')) {
                 return name
               }
               // Skip everything else (raw input FASTQ files)
               return null
             }

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(reads),
          val(condition), val(timepoint), val(replicate)
    val   data_type

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path('final_R1.fastq'),
          path('final_R2.fastq'),
          val(condition), val(timepoint), val(replicate),
          emit: cleaned_reads

    path 'fastqc_raw/*',                  emit: fastqc_raw
    path 'fastqc_final/*',                emit: fastqc_final
    path 'trim_stats.tsv',                emit: trim_stats
    path 'umi_stats.tsv',    optional: true, emit: umi_stats
    path 'barcode_umi_profile.json', optional: true, emit: barcode_umi_profile
    path 'barcode_umi_detect.log',   optional: true, emit: barcode_umi_detect_log
    path 'cutadapt_report.txt',           emit: cutadapt_report
    path 'umi_extract.log',  optional: true, emit: umi_log
    path 'preprocess_reads.log',          emit: log
    path 'README_01_trimmed_fastq.txt'

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  # umi_tools imports matplotlib; avoid font-cache stall
  export MPLCONFIGDIR="\${TMPDIR:-/tmp}/matplotlib"

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  # tee stdout (not plain redirect) so .command.out also receives the live
  # "PREP |" progress — otherwise the monitor has nothing to tail for this process.
  exec > >(tee preprocess_reads.log)
  exec 2> >(tee -a preprocess_reads.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "preprocess_and_quality_filter_reads" "Unexpected process failure" "Check preprocess_reads.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "PREP | START | sample=${sample_id} | mode=${data_type ?: "SE"} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID='${sample_id}'
  THREADS=${task.cpus}
  MODE='${(data_type ?: "SE").toString()}'
  
  R1='${reads[0]}'
  R2='${(reads.size() > 1 ? reads[1] : "")}'

  echo "PREP | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "PREP | CONFIG | Mode: \${MODE}"
  echo "PREP | CONFIG | Threads: \${THREADS}"
  echo "PREP | CONFIG | R1: \${R1}"
  if [[ -n "\${R2}" ]]; then
    echo "PREP | CONFIG | R2: \${R2}"
  fi

  # ── Adapter Trimming Parameters ──
  TRIM_ENABLED=\$([ '${params.adapter_trimming?.enabled == true}' = 'true' ] && echo 1 || echo 0)
  ADAPTER1='${(params.adapter_trimming?.adapter1 ?: "").toString()}'
  ADAPTER2='${(params.adapter_trimming?.adapter2 ?: "").toString()}'
  USER_MINLEN=${(params.adapter_trimming?.minlen ?: 0)}

  echo "PREP | CONFIG | Adapter trimming: \$([ \${TRIM_ENABLED} -eq 1 ] && echo "enabled" || echo "disabled")"
  if [[ \${TRIM_ENABLED} -eq 1 ]]; then
    echo "PREP | CONFIG | Adapter 1: \${ADAPTER1:-none}"
    echo "PREP | CONFIG | Adapter 2: \${ADAPTER2:-none}"
  fi

  # ── Barcode Parameters (primary + optional second) ──
  # Each barcode slot is independent: which READ it's on (R1/R2) x which
  # END (5'/3') x LENGTH x how much to trust the stated location (detect_mode).
  BC1_ENABLED=\$([ '${params.barcode?.enabled == true}' = 'true' ] && echo 1 || echo 0)
  BC1_LENGTH=${(params.barcode?.length ?: 0)}
  BC1_LOCATION='${(params.barcode?.location ?: "5").toString()}'
  BC1_READ='${(params.barcode?.read ?: "R1").toString().toUpperCase()}'
  BC1_DETECT_MODE='${(params.barcode?.detect_mode ?: "trust").toString().toLowerCase()}'
  BC1_DETECT_TOL=${(params.barcode?.detect_tolerance ?: 1)}

  BC2_ENABLED=\$([ '${params.barcode?.enabled2 == true}' = 'true' ] && echo 1 || echo 0)
  BC2_LENGTH=${(params.barcode?.length2 ?: 0)}
  BC2_LOCATION='${(params.barcode?.location2 ?: "5").toString()}'
  BC2_READ='${(params.barcode?.read2 ?: "R2").toString().toUpperCase()}'
  BC2_DETECT_MODE='${(params.barcode?.detect_mode2 ?: "trust").toString().toLowerCase()}'
  BC2_DETECT_TOL=${(params.barcode?.detect_tolerance2 ?: 1)}

  echo "PREP | CONFIG | Barcode 1: \$([ \${BC1_ENABLED} -eq 1 ] && echo "enabled (\${BC1_READ}, \${BC1_LOCATION}', \${BC1_LENGTH}bp, mode=\${BC1_DETECT_MODE})" || echo "disabled")"
  if [[ \${BC2_ENABLED} -eq 1 ]]; then
    echo "PREP | CONFIG | Barcode 2: enabled (\${BC2_READ}, \${BC2_LOCATION}', \${BC2_LENGTH}bp, mode=\${BC2_DETECT_MODE})"
  fi

  # ── UMI Parameters ──
  UMI_ENABLED=\$([ '${params.umi?.enabled == true}' = 'true' ] && echo 1 || echo 0)
  UMI_LENGTH=${(params.umi?.length ?: 0)}
  UMI_LOCATION='${(params.umi?.location ?: "5").toString()}'
  UMI_READ='${(params.umi?.read ?: "R1").toString().toUpperCase()}'
  UMI_DETECT_MODE='${(params.umi?.detect_mode ?: "trust").toString().toLowerCase()}'
  UMI_DETECT_TOL=${(params.umi?.detect_tolerance ?: 1)}

  echo "PREP | CONFIG | UMI: \$([ \${UMI_ENABLED} -eq 1 ] && echo "enabled (\${UMI_READ}, \${UMI_LOCATION}', \${UMI_LENGTH}bp, mode=\${UMI_DETECT_MODE})" || echo "disabled")"

  if [[ "\${MODE}" != "PE" ]]; then
    # No R2 exists in SE mode: a slot pointed at R2 is disabled (with a
    # warning) rather than silently reinterpreted as R1 -- guessing which
    # end the user actually meant would risk trimming the wrong bases.
    if [[ "\${BC1_READ}" == "R2" && \${BC1_ENABLED} -eq 1 ]]; then
      echo "PREP | CONFIG | WARNING | Barcode 1 targets R2 but mode is SE (no R2 exists) -- disabling"
      BC1_ENABLED=0
    fi
    if [[ "\${BC2_READ}" == "R2" && \${BC2_ENABLED} -eq 1 ]]; then
      echo "PREP | CONFIG | WARNING | Barcode 2 targets R2 but mode is SE (no R2 exists) -- disabling"
      BC2_ENABLED=0
    fi
    if [[ "\${UMI_READ}" == "R2" && \${UMI_ENABLED} -eq 1 ]]; then
      echo "PREP | CONFIG | WARNING | UMI targets R2 but mode is SE (no R2 exists) -- disabling"
      UMI_ENABLED=0
    fi
  fi

  # ── QC Parameters ──
  QC_ENABLED=\$([ '${(params.qc?.enabled == null ? true : params.qc?.enabled)}' = 'true' ] && echo 1 || echo 0)
  QC_RAW=\$([ '${(params.fastqc_raw == null ? true : params.fastqc_raw)}' = 'true' ] && echo 1 || echo 0)
  FINAL_MINLEN=${(params.advanced?.prep_insert_minlen ?: 12)}

  echo "PREP | CONFIG | QC (FastQC): \$([ \${QC_ENABLED} -eq 1 ] && echo "enabled" || echo "disabled")"
  echo "PREP | CONFIG | QC on raw reads: \$([ \${QC_RAW} -eq 1 ] && echo "yes" || echo "no")"
  echo "PREP | CONFIG | Final minimum length: \${FINAL_MINLEN}bp"

  ###########################################################################
  # 1b) BARCODE/UMI QC DETECTION (only runs if any slot asks for it)
  ###########################################################################
  # "trust" (default) never reaches the detector call below -- zero cost,
  # identical behavior to before this feature existed. "auto"/"verify" on
  # ANY slot triggers ONE scan (bin/detect_barcode_umi.py) of a subsample of
  # the RAW reads, comparing the letter-mix at each end to the letter-mix in
  # the middle of the read; every slot's decision below reuses that single
  # scan via --profile-in (near-instant), so multiple auto/verify slots
  # don't each pay for a fresh FASTQ pass. See the script's own docstring
  # for the full method and its documented limitations (short tags,
  # near-uniform-GC organisms).

  NEEDS_DETECT=0
  for m in "\${BC1_DETECT_MODE}" "\${BC2_DETECT_MODE}" "\${UMI_DETECT_MODE}"; do
    [[ "\$m" == "auto" || "\$m" == "verify" ]] && NEEDS_DETECT=1
  done

  if [[ \${NEEDS_DETECT} -eq 1 ]]; then
    echo "PREP | DETECT | At least one barcode/UMI slot uses auto or verify -- running QC scan..."
    : > barcode_umi_detect.log
    python3 detect_barcode_umi.py --r1 "\${R1}" \${R2:+--r2 "\${R2}"} \\
      --n-reads 200000 --out barcode_umi_profile.json --emit-shell \\
      > detect_vars.sh 2>>barcode_umi_detect.log
    cat detect_vars.sh >> barcode_umi_detect.log
    source detect_vars.sh   # -> BEST_BARCODE_READ/LOCATION/LENGTH, BEST_UMI_READ/LOCATION/LENGTH

    # ── Barcode 1: auto ──
    if [[ "\${BC1_DETECT_MODE}" == "auto" && \${BC1_ENABLED} -eq 1 ]]; then
      if [[ -z "\${BEST_BARCODE_READ}" || "\${BEST_BARCODE_LENGTH:-0}" -le 0 ]]; then
        echo "PREP | DETECT | Barcode 1: mode=auto, QC found nothing conclusive -> disabling" | tee -a barcode_umi_detect.log
        BC1_ENABLED=0
      else
        echo "PREP | DETECT | Barcode 1: mode=auto, QC proposes \${BEST_BARCODE_READ} \${BEST_BARCODE_LOCATION}' \${BEST_BARCODE_LENGTH}bp (was: \${BC1_READ} \${BC1_LOCATION}' \${BC1_LENGTH}bp)" | tee -a barcode_umi_detect.log
        BC1_READ="\${BEST_BARCODE_READ}"; BC1_LOCATION="\${BEST_BARCODE_LOCATION}"; BC1_LENGTH="\${BEST_BARCODE_LENGTH}"
      fi
    fi

    # ── Barcode 1: verify ──
    if [[ "\${BC1_DETECT_MODE}" == "verify" && \${BC1_ENABLED} -eq 1 ]]; then
      # NOTE: the check runs as the condition of this `if` (rather than
      # `VAR=$(...); RC=$?` as two separate statements) on purpose -- under
      # `set -e`, a plain assignment statement that fails aborts the script
      # right there via the ERR trap, before the next line even runs, which
      # would skip the specific tracktx_error message below in favor of a
      # generic "Unexpected process failure". Bash exempts a command's exit
      # status from triggering errexit when that command IS the `if` test.
      if BC1_CHECK=\$(python3 detect_barcode_umi.py --profile-in barcode_umi_profile.json \\
        --check-kind barcode --check-read "\${BC1_READ}" --check-location "\${BC1_LOCATION}" \\
        --check-length "\${BC1_LENGTH}" --tolerance "\${BC1_DETECT_TOL}" 2>>barcode_umi_detect.log); then
        BC1_CHECK_RC=0
      else
        BC1_CHECK_RC=\$?
      fi
      echo "\${BC1_CHECK}" >> barcode_umi_detect.log
      if [[ \${BC1_CHECK_RC} -ne 0 ]]; then
        BC1_ELSEWHERE=\$(python3 - "\${BC1_CHECK}" <<'PYEOF' 2>/dev/null || echo "unknown"
import json, sys
d = json.loads(sys.argv[1])
b = d.get("best_guess_elsewhere")
if b:
    print(f"{b['read']} {b['location']}' {b['length']}bp")
else:
    print("nothing conclusive found elsewhere either")
PYEOF
)
        tracktx_error "preprocess_and_quality_filter_reads" \\
          "Barcode 1 mismatch: you said \${BC1_READ} \${BC1_LOCATION}' \${BC1_LENGTH}bp, but the QC scan of the raw reads doesn't support that (nearest match: \${BC1_ELSEWHERE})" \\
          "Fix params.barcode.read/location/length, widen params.barcode.detect_tolerance, or switch params.barcode.detect_mode to 'auto'. Full profile: barcode_umi_profile.json"
      else
        echo "PREP | DETECT | Barcode 1: mode=verify, QC confirms \${BC1_READ} \${BC1_LOCATION}' \${BC1_LENGTH}bp" | tee -a barcode_umi_detect.log
      fi
    fi

    # ── Barcode 2: auto (excludes barcode 1's corner so it isn't proposed twice) ──
    if [[ "\${BC2_DETECT_MODE}" == "auto" && \${BC2_ENABLED} -eq 1 ]]; then
      python3 detect_barcode_umi.py --profile-in barcode_umi_profile.json \\
        --exclude-read "\${BC1_READ}" --exclude-location "\${BC1_LOCATION}" \\
        --emit-shell --out /dev/null > detect_vars_bc2.sh 2>>barcode_umi_detect.log
      cat detect_vars_bc2.sh >> barcode_umi_detect.log
      BC2_BEST_READ=\$(grep '^BEST_BARCODE_READ=' detect_vars_bc2.sh | cut -d= -f2)
      BC2_BEST_LOC=\$(grep '^BEST_BARCODE_LOCATION=' detect_vars_bc2.sh | cut -d= -f2)
      BC2_BEST_LEN=\$(grep '^BEST_BARCODE_LENGTH=' detect_vars_bc2.sh | cut -d= -f2)
      if [[ -z "\${BC2_BEST_READ}" || "\${BC2_BEST_LEN:-0}" -le 0 ]]; then
        echo "PREP | DETECT | Barcode 2: mode=auto, QC found nothing conclusive (beyond barcode 1's corner) -> disabling" | tee -a barcode_umi_detect.log
        BC2_ENABLED=0
      else
        echo "PREP | DETECT | Barcode 2: mode=auto, QC proposes \${BC2_BEST_READ} \${BC2_BEST_LOC}' \${BC2_BEST_LEN}bp (was: \${BC2_READ} \${BC2_LOCATION}' \${BC2_LENGTH}bp)" | tee -a barcode_umi_detect.log
        BC2_READ="\${BC2_BEST_READ}"; BC2_LOCATION="\${BC2_BEST_LOC}"; BC2_LENGTH="\${BC2_BEST_LEN}"
      fi
    fi

    # ── Barcode 2: verify ──
    if [[ "\${BC2_DETECT_MODE}" == "verify" && \${BC2_ENABLED} -eq 1 ]]; then
      if BC2_CHECK=\$(python3 detect_barcode_umi.py --profile-in barcode_umi_profile.json \\
        --check-kind barcode --check-read "\${BC2_READ}" --check-location "\${BC2_LOCATION}" \\
        --check-length "\${BC2_LENGTH}" --tolerance "\${BC2_DETECT_TOL}" 2>>barcode_umi_detect.log); then
        BC2_CHECK_RC=0
      else
        BC2_CHECK_RC=\$?
      fi
      echo "\${BC2_CHECK}" >> barcode_umi_detect.log
      if [[ \${BC2_CHECK_RC} -ne 0 ]]; then
        BC2_ELSEWHERE=\$(python3 - "\${BC2_CHECK}" <<'PYEOF' 2>/dev/null || echo "unknown"
import json, sys
d = json.loads(sys.argv[1])
b = d.get("best_guess_elsewhere")
if b:
    print(f"{b['read']} {b['location']}' {b['length']}bp")
else:
    print("nothing conclusive found elsewhere either")
PYEOF
)
        tracktx_error "preprocess_and_quality_filter_reads" \\
          "Barcode 2 mismatch: you said \${BC2_READ} \${BC2_LOCATION}' \${BC2_LENGTH}bp, but the QC scan of the raw reads doesn't support that (nearest match: \${BC2_ELSEWHERE})" \\
          "Fix params.barcode.read2/location2/length2, widen params.barcode.detect_tolerance2, or switch params.barcode.detect_mode2 to 'auto'. Full profile: barcode_umi_profile.json"
      else
        echo "PREP | DETECT | Barcode 2: mode=verify, QC confirms \${BC2_READ} \${BC2_LOCATION}' \${BC2_LENGTH}bp" | tee -a barcode_umi_detect.log
      fi
    fi

    # ── UMI: auto ──
    if [[ "\${UMI_DETECT_MODE}" == "auto" && \${UMI_ENABLED} -eq 1 ]]; then
      if [[ -z "\${BEST_UMI_READ}" || "\${BEST_UMI_LENGTH:-0}" -le 0 ]]; then
        echo "PREP | DETECT | UMI: mode=auto, QC found nothing conclusive -> disabling" | tee -a barcode_umi_detect.log
        UMI_ENABLED=0
      else
        echo "PREP | DETECT | UMI: mode=auto, QC proposes \${BEST_UMI_READ} \${BEST_UMI_LOCATION}' \${BEST_UMI_LENGTH}bp (was: \${UMI_READ} \${UMI_LOCATION}' \${UMI_LENGTH}bp)" | tee -a barcode_umi_detect.log
        UMI_READ="\${BEST_UMI_READ}"; UMI_LOCATION="\${BEST_UMI_LOCATION}"; UMI_LENGTH="\${BEST_UMI_LENGTH}"
      fi
    fi

    # ── UMI: verify ──
    if [[ "\${UMI_DETECT_MODE}" == "verify" && \${UMI_ENABLED} -eq 1 ]]; then
      if UMI_CHECK=\$(python3 detect_barcode_umi.py --profile-in barcode_umi_profile.json \\
        --check-kind umi --check-read "\${UMI_READ}" --check-location "\${UMI_LOCATION}" \\
        --check-length "\${UMI_LENGTH}" --tolerance "\${UMI_DETECT_TOL}" 2>>barcode_umi_detect.log); then
        UMI_CHECK_RC=0
      else
        UMI_CHECK_RC=\$?
      fi
      echo "\${UMI_CHECK}" >> barcode_umi_detect.log
      if [[ \${UMI_CHECK_RC} -ne 0 ]]; then
        UMI_ELSEWHERE=\$(python3 - "\${UMI_CHECK}" <<'PYEOF' 2>/dev/null || echo "unknown"
import json, sys
d = json.loads(sys.argv[1])
b = d.get("best_guess_elsewhere")
if b:
    print(f"{b['read']} {b['location']}' {b['length']}bp")
else:
    print("nothing conclusive found elsewhere either")
PYEOF
)
        tracktx_error "preprocess_and_quality_filter_reads" \\
          "UMI mismatch: you said \${UMI_READ} \${UMI_LOCATION}' \${UMI_LENGTH}bp, but the QC scan of the raw reads doesn't support that (nearest match: \${UMI_ELSEWHERE})" \\
          "Fix params.umi.read/location/length, widen params.umi.detect_tolerance, or switch params.umi.detect_mode to 'auto'. Full profile: barcode_umi_profile.json"
      else
        echo "PREP | DETECT | UMI: mode=verify, QC confirms \${UMI_READ} \${UMI_LOCATION}' \${UMI_LENGTH}bp" | tee -a barcode_umi_detect.log
      fi
    fi

  else
    echo "PREP | DETECT | All slots on detect_mode=trust -- skipping QC scan (single-pass default)"
  fi

  # ── Calculate Length Requirements ──
  # Barcode loss is tallied PER READ (each slot can independently sit on R1
  # or R2 now) for logging only -- it does NOT feed the cutadapt -m floor
  # below, because cutadapt applies -m to each read's length AFTER its own
  # -u/-U cuts in the same invocation, so barcode bases are already gone by
  # the time -m is checked; adding BARCODE_LOSS to the floor would double
  # count it. (Previously R1's and R2's barcode lengths were simply summed
  # together into one number regardless of which read they were actually on
  # -- harmless when barcode1 was always R1 and barcode2 was always R2, but
  # not generally correct once either slot can point at either read.)
  R1_BC_LOSS=0; R2_BC_LOSS=0
  [[ \${BC1_ENABLED} -eq 1 ]] && { [[ "\${BC1_READ}" == "R2" ]] && R2_BC_LOSS=\$((R2_BC_LOSS + BC1_LENGTH)) || R1_BC_LOSS=\$((R1_BC_LOSS + BC1_LENGTH)); }
  [[ \${BC2_ENABLED} -eq 1 ]] && { [[ "\${BC2_READ}" == "R2" ]] && R2_BC_LOSS=\$((R2_BC_LOSS + BC2_LENGTH)) || R1_BC_LOSS=\$((R1_BC_LOSS + BC2_LENGTH)); }

  # params.adapter_trimming.minlen ("user override") used to be computed
  # into EFFECTIVE_MINLEN and then never applied anywhere -- the cutadapt
  # -m floor actually used was FINAL_MINLEN (params.advanced.prep_insert_minlen,
  # default 12) alone, so a params.yaml setting minlen above/below 12 had no
  # effect. Fixed: the user override now raises (never lowers) the floor,
  # same UMI_LOSS handling as before (UMI hasn't been extracted yet at this
  # point -- that happens later, via umi_tools -- so its bases are still in
  # the read when cutadapt's -m check runs).
  UMI_LOSS=\$(( UMI_ENABLED == 1 ? UMI_LENGTH : 0 ))
  BASE_MINLEN=\$(( USER_MINLEN > FINAL_MINLEN ? USER_MINLEN : FINAL_MINLEN ))
  PRE_UMI_MINLEN=\$(( BASE_MINLEN + UMI_LOSS ))

  echo "PREP | CONFIG | Length calculations:"
  echo "PREP | CONFIG |   Barcode loss (R1): \${R1_BC_LOSS}bp, (R2): \${R2_BC_LOSS}bp (informational -- already removed by cutadapt before -m applies)"
  echo "PREP | CONFIG |   UMI loss: \${UMI_LOSS}bp (on \${UMI_READ})"
  echo "PREP | CONFIG |   Default min length: \${FINAL_MINLEN}bp"
  echo "PREP | CONFIG |   User min length override: \${USER_MINLEN}bp"
  echo "PREP | CONFIG |   Base min length applied: \${BASE_MINLEN}bp"
  echo "PREP | CONFIG |   Min length filter (cutadapt -m, pre-UMI-extraction): \${PRE_UMI_MINLEN}bp"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "PREP | VALIDATE | Checking inputs..."

  # Fix: Files may be gzipped but have .fastq extension.
  # Case 1: Symlink from download_sra_samples (SRR_R1.fastq -> SRR_R1.fastq.gz). Use target, do NOT mv
  #   (mv would overwrite the real .gz file with the symlink and break it).
  # Case 2: Regular file with .fastq extension but gzip magic bytes. Rename to .gz.
  # IMPORTANT: This function is used with \$(fix_gzip_extension "\${R1}"). Only the final
  # path must go to stdout; all log messages must use >&2.
  fix_gzip_extension() {
    local f="\$1"
    [[ -z "\$f" || ! -e "\$f" ]] && echo "\$f" && return
    [[ "\$f" == *.gz ]] && echo "\$f" && return
    # Symlink: use target (the actual .gz file) - do not mv
    if [[ -L "\$f" ]]; then
      local dir base target
      dir=\$(dirname "\$f")
      base=\$(basename "\$f")
      target=\$(readlink "\$f")  # target may be relative
      if [[ "\$target" != /* ]]; then
        target="\${dir}/\${target}"
      fi
      if [[ -f "\$target" && "\$target" == *.gz ]]; then
        echo "PREP | VALIDATE | Detected symlink .fastq -> .gz, using target: \$target" >&2
        echo "\$target"
        return
      fi
    fi
    # Regular file: check gzip magic bytes
    local magic
    magic=\$(head -c 2 "\$f" 2>/dev/null | od -A n -t x1 2>/dev/null | tr -d ' \\n' | head -c 4)
    if [[ "\$magic" == "1f8b" ]]; then
      echo "PREP | VALIDATE | Detected gzipped content with .fastq extension, renaming to .gz: \$f" >&2
      mv "\$f" "\${f}.gz"
      echo "\${f}.gz"
    else
      echo "\$f"
    fi
  }
  # fix_gzip_extension returns path via stdout; use tail -1 so only the path is used
  # (defensive: if any log line were echoed to stdout, the path is always last)
  R1_FIXED=\$(fix_gzip_extension "\${R1}")
  R1=\$(echo "\${R1_FIXED}" | tail -1)
  if [[ -n "\${R2}" && -f "\${R2}" ]]; then
    R2_FIXED=\$(fix_gzip_extension "\${R2}")
    R2=\$(echo "\${R2_FIXED}" | tail -1)
  fi

  if [[ "\${MODE}" != "SE" && "\${MODE}" != "PE" ]]; then
    tracktx_error "preprocess_and_quality_filter_reads" "Mode must be SE or PE, got: \${MODE}" "Check data_type parameter"
  fi

  if [[ "\${MODE}" == "PE" && -z "\${R2}" ]]; then
    tracktx_error "preprocess_and_quality_filter_reads" "Paired-end mode but R2 file is missing" "Add file2 to samplesheet for PE samples"
  fi

  if [[ ! -f "\${R1}" ]]; then
    tracktx_error "preprocess_and_quality_filter_reads" "R1 file not found: \${R1}" "Check samplesheet file1 paths"
  fi

  if [[ "\${MODE}" == "PE" && ! -f "\${R2}" ]]; then
    tracktx_error "preprocess_and_quality_filter_reads" "R2 file not found: \${R2}" "Check samplesheet file2 paths"
  fi

  # Check file sizes
  R1_SIZE=\$(stat -c%s "\${R1}" 2>/dev/null || stat -f%z "\${R1}" 2>/dev/null || echo "unknown")
  echo "PREP | VALIDATE | R1 size: \${R1_SIZE} bytes"
  
  if [[ "\${MODE}" == "PE" ]]; then
    R2_SIZE=\$(stat -c%s "\${R2}" 2>/dev/null || stat -f%z "\${R2}" 2>/dev/null || echo "unknown")
    echo "PREP | VALIDATE | R2 size: \${R2_SIZE} bytes"
  fi

  echo "PREP | VALIDATE | Input validation complete"

  ###########################################################################
  # 3) HELPER FUNCTIONS
  ###########################################################################

  # Count reads in FASTQ file (handles .gz)
  count_reads() {
    local file="\$1"
    if [[ "\$file" == *.gz ]]; then
      gzip -cd "\$file" | awk 'END{print NR/4}' 2>/dev/null || echo 0
    else
      awk 'END{print NR/4}' "\$file" 2>/dev/null || echo 0
    fi
  }

  ###########################################################################
  # 4) QC ON RAW READS (Optional)
  ###########################################################################

  mkdir -p fastqc_raw

  if [[ \${QC_RAW} -eq 1 && \${QC_ENABLED} -eq 1 ]]; then
    echo "PREP | QC-RAW | Running FastQC on raw reads..."
    
    if [[ "\${MODE}" == "PE" ]]; then
      fastqc --quiet --threads "\${THREADS}" -o fastqc_raw "\${R1}" "\${R2}"
    else
      fastqc --quiet --threads "\${THREADS}" -o fastqc_raw "\${R1}"
    fi
    
    echo "PREP | QC-RAW | FastQC reports generated"
  else
    echo "PREP | QC-RAW | Skipped (disabled)"
    echo "QC on raw reads disabled." > fastqc_raw/README.txt
  fi

  ###########################################################################
  # 5) SINGLE-PASS CUTADAPT
  ###########################################################################

  echo "PREP | CUTADAPT | Running single-pass trimming..."
  echo "PREP | CUTADAPT | Output: preumi_R*.fastq"

  # Build cutadapt command
  CUTADAPT_CMD=("-j" "\${THREADS}")
  
  # Minimum length filter (pre-UMI)
  CUTADAPT_CMD+=("-m" "\${PRE_UMI_MINLEN}")
  echo "PREP | CUTADAPT | Minimum length filter: \${PRE_UMI_MINLEN}bp"

  # Adapter trimming
  if [[ \${TRIM_ENABLED} -eq 1 ]]; then
    if [[ -n "\${ADAPTER1}" ]]; then
      CUTADAPT_CMD+=("-a" "\${ADAPTER1}")
      echo "PREP | CUTADAPT | R1 adapter: \${ADAPTER1}"
    fi
    if [[ "\${MODE}" == "PE" && -n "\${ADAPTER2}" ]]; then
      CUTADAPT_CMD+=("-A" "\${ADAPTER2}")
      echo "PREP | CUTADAPT | R2 adapter: \${ADAPTER2}"
    fi
  fi

  # Barcode removal — each slot is routed to cutadapt's R1 (-u) or R2 (-U)
  # flag by its OWN \${BCn_READ} setting, not by slot number. cutadapt allows
  # -u (or -U) to appear twice for the same read -- once positive (5' cut),
  # once negative (3' cut) -- so two slots landing on the same read at
  # opposite ends both apply correctly in one invocation.
  if [[ \${BC1_ENABLED} -eq 1 && \${BC1_LENGTH} -gt 0 ]]; then
    if [[ "\${BC1_READ}" == "R2" ]]; then
      [[ "\${BC1_LOCATION}" == "3" ]] && CUTADAPT_CMD+=("-U" "-\${BC1_LENGTH}") || CUTADAPT_CMD+=("-U" "\${BC1_LENGTH}")
      echo "PREP | CUTADAPT | Barcode 1: \${BC1_LENGTH}bp from R2 \${BC1_LOCATION}' end"
    else
      [[ "\${BC1_LOCATION}" == "3" ]] && CUTADAPT_CMD+=("-u" "-\${BC1_LENGTH}") || CUTADAPT_CMD+=("-u" "\${BC1_LENGTH}")
      echo "PREP | CUTADAPT | Barcode 1: \${BC1_LENGTH}bp from R1 \${BC1_LOCATION}' end"
    fi
  fi

  if [[ \${BC2_ENABLED} -eq 1 && \${BC2_LENGTH} -gt 0 ]]; then
    if [[ "\${BC2_READ}" == "R2" ]]; then
      [[ "\${BC2_LOCATION}" == "3" ]] && CUTADAPT_CMD+=("-U" "-\${BC2_LENGTH}") || CUTADAPT_CMD+=("-U" "\${BC2_LENGTH}")
      echo "PREP | CUTADAPT | Barcode 2: \${BC2_LENGTH}bp from R2 \${BC2_LOCATION}' end"
    else
      [[ "\${BC2_LOCATION}" == "3" ]] && CUTADAPT_CMD+=("-u" "-\${BC2_LENGTH}") || CUTADAPT_CMD+=("-u" "\${BC2_LENGTH}")
      echo "PREP | CUTADAPT | Barcode 2: \${BC2_LENGTH}bp from R1 \${BC2_LOCATION}' end"
    fi
  fi

  # Run cutadapt
  echo "PREP | CUTADAPT | Processing reads..."
  if [[ "\${MODE}" == "PE" ]]; then
    cutadapt "\${CUTADAPT_CMD[@]}" \\
             -o preumi_R1.fastq \\
             -p preumi_R2.fastq \\
             "\${R1}" "\${R2}" | tee cutadapt_report.txt
  else
    cutadapt "\${CUTADAPT_CMD[@]}" \\
             -o preumi_R1.fastq \\
             "\${R1}" | tee cutadapt_report.txt
  fi

  echo "PREP | CUTADAPT | Trimming complete"

  ###########################################################################
  # 6) UMI EXTRACTION
  ###########################################################################

  if [[ \${UMI_ENABLED} -eq 1 && \${UMI_LENGTH} -gt 0 ]]; then
    echo "PREP | UMI | Extracting UMI sequences..."
    echo "PREP | UMI | Read: \${UMI_READ}"
    echo "PREP | UMI | Length: \${UMI_LENGTH}bp"
    echo "PREP | UMI | Location: \${UMI_LOCATION}' end"

    # Build UMI pattern (e.g., NNNNNN for 6bp UMI)
    UMI_PATTERN=\$(printf 'N%.0s' \$(seq 1 \${UMI_LENGTH}))

    # Location flag (umi_tools applies this to whichever pattern is given below)
    UMI_END_FLAG=""
    if [[ "\${UMI_LOCATION}" == "3" ]]; then
      UMI_END_FLAG="--3prime"
    fi

    # Which read the UMI lives on decides --bc-pattern (R1) vs --bc-pattern2
    # (R2). umi_tools accepts EITHER on its own for paired-end input -- it
    # does not require a pattern on both reads (confirmed against umi_tools
    # 1.1.5's own option validation: "Must supply --bc-pattern and/or
    # --bc-pattern2 if paired-end"). This is what unblocks a library like
    # the E. coli MG1655 cohort entry, whose UMI is on R2 while every prior
    # cohort entry only ever had one on R1.
    UMI_PATTERN_FLAG="--bc-pattern=\${UMI_PATTERN}"
    if [[ "\${MODE}" == "PE" && "\${UMI_READ}" == "R2" ]]; then
      UMI_PATTERN_FLAG="--bc-pattern2=\${UMI_PATTERN}"
    fi

    # Extract UMI, parallelized across read chunks: umi_tools extract is a
    # stateless per-read regex operation with no cross-read state, and has
    # no thread flag of its own (same limitation as umi_tools dedup, fixed
    # the same way in module 06). Splitting preumi_R{1,2}.fastq into N
    # equal, record-aligned chunks (same line-count boundaries for both
    # mates, so R1/R2 pairing stays intact), running extract on each chunk
    # independently, then concatenating outputs back in the same chunk
    # order produces output identical to one serial call.
    mkdir -p umi_chunks
    rm -f umi_chunks/* umi_extract_fail.flag 2>/dev/null || true

    R1_LINES=\$(wc -l < preumi_R1.fastq | tr -d ' ')
    CHUNK_LINES=\$(( (R1_LINES / THREADS / 4 + 1) * 4 ))
    [ "\${CHUNK_LINES}" -lt 4 ] && CHUNK_LINES=4

    split -d -a 4 -l "\${CHUNK_LINES}" preumi_R1.fastq umi_chunks/r1_
    if [[ "\${MODE}" == "PE" ]]; then
      split -d -a 4 -l "\${CHUNK_LINES}" preumi_R2.fastq umi_chunks/r2_
    fi

    extract_chunk() {
      local r1_in="\$1" r1_out="\$2" r2_in="\$3" r2_out="\$4"
      if [[ -n "\${r2_in}" ]]; then
        if ! umi_tools extract \${UMI_PATTERN_FLAG} \${UMI_END_FLAG} \\
             -I "\${r1_in}" -S "\${r1_out}" \\
             --read2-in "\${r2_in}" --read2-out "\${r2_out}" \\
             --log="\${r1_out}.log" >>"\${r1_out}.err" 2>&1; then
          echo "FAIL: \${r1_in}" >> umi_extract_fail.flag
        fi
      else
        if ! umi_tools extract \${UMI_PATTERN_FLAG} \${UMI_END_FLAG} \\
             -I "\${r1_in}" -S "\${r1_out}" \\
             --log="\${r1_out}.log" >>"\${r1_out}.err" 2>&1; then
          echo "FAIL: \${r1_in}" >> umi_extract_fail.flag
        fi
      fi
    }

    EXTRACT_PIDS=()
    launch_extract() {
      while :; do
        local alive=0 p
        for p in "\${EXTRACT_PIDS[@]:-}"; do
          [ -n "\${p}" ] && kill -0 "\${p}" 2>/dev/null && alive=\$(( alive + 1 ))
        done
        [ "\${alive}" -lt "\${THREADS}" ] && break
        sleep 0.5
      done
      extract_chunk "\$1" "\$2" "\$3" "\$4" &
      EXTRACT_PIDS+=(\$!)
    }

    R1_OUT_CHUNKS=()
    R2_OUT_CHUNKS=()
    for r1chunk in umi_chunks/r1_*; do
      suffix="\${r1chunk#umi_chunks/r1_}"
      r1out="umi_chunks/out_r1_\${suffix}.fastq"
      if [[ "\${MODE}" == "PE" ]]; then
        r2chunk="umi_chunks/r2_\${suffix}"
        r2out="umi_chunks/out_r2_\${suffix}.fastq"
        launch_extract "\${r1chunk}" "\${r1out}" "\${r2chunk}" "\${r2out}"
        R2_OUT_CHUNKS+=("\${r2out}")
      else
        launch_extract "\${r1chunk}" "\${r1out}" "" ""
      fi
      R1_OUT_CHUNKS+=("\${r1out}")
    done

    for p in "\${EXTRACT_PIDS[@]}"; do
      wait "\${p}" || true
    done

    if [[ -s umi_extract_fail.flag ]]; then
      echo "PREP | ERROR | Failed UMI extract chunk(s):"
      sed 's/^/PREP | ERROR |   /' umi_extract_fail.flag
      cat umi_chunks/*.err >&2 2>/dev/null || true
      tracktx_error "preprocess_and_quality_filter_reads" "umi_tools extraction failed" "Check preprocess_reads.log in work dir"
    fi

    cat "\${R1_OUT_CHUNKS[@]}" > final_R1.fastq
    if [[ "\${MODE}" == "PE" ]]; then
      cat "\${R2_OUT_CHUNKS[@]}" > final_R2.fastq
    else
      # Create empty R2 for SE (consistent tuple shape)
      : > final_R2.fastq
    fi
    cat umi_chunks/*.log > umi_extract.log 2>/dev/null || true
    rm -rf umi_chunks

    echo "PREP | UMI | Extraction complete"

  else
    echo "PREP | UMI | UMI extraction disabled, passing through..."
    
    if [[ "\${MODE}" == "PE" ]]; then
      cp -f preumi_R1.fastq final_R1.fastq
      cp -f preumi_R2.fastq final_R2.fastq
    else
      cp -f preumi_R1.fastq final_R1.fastq
      : > final_R2.fastq  # Empty stub for SE
    fi

    echo "PREP | UMI | Pass-through complete"
  fi

  ###########################################################################
  # 7) QC ON FINAL READS
  ###########################################################################

  mkdir -p fastqc_final

  if [[ \${QC_ENABLED} -eq 1 ]]; then
    echo "PREP | QC-FINAL | Running FastQC on final reads..."
    
    if [[ "\${MODE}" == "PE" ]]; then
      fastqc --quiet --threads "\${THREADS}" -o fastqc_final final_R1.fastq final_R2.fastq
    else
      fastqc --quiet --threads "\${THREADS}" -o fastqc_final final_R1.fastq
    fi
    
    echo "PREP | QC-FINAL | FastQC reports generated"
  else
    echo "PREP | QC-FINAL | Skipped (QC disabled)"
    echo "QC disabled." > fastqc_final/README.txt
  fi

  ###########################################################################
  # 8) COMPUTE STATISTICS
  ###########################################################################

  echo "PREP | STATS | Computing read statistics..."

  # Count reads at each stage
  R1_IN=\$(count_reads "\${R1}")
  R1_PRE=\$(count_reads preumi_R1.fastq)
  R1_FINAL=\$(count_reads final_R1.fastq)

  if [[ "\${MODE}" == "PE" ]]; then
    R2_IN=\$(count_reads "\${R2}")
    R2_PRE=\$(count_reads preumi_R2.fastq)
    R2_FINAL=\$(count_reads final_R2.fastq)
  else
    R2_IN=0
    R2_PRE=0
    R2_FINAL=0
  fi

  echo "PREP | STATS | R1: \${R1_IN} → \${R1_PRE} → \${R1_FINAL} reads"
  if [[ "\${MODE}" == "PE" ]]; then
    echo "PREP | STATS | R2: \${R2_IN} → \${R2_PRE} → \${R2_FINAL} reads"
  fi

  # Trimming statistics
  {
    echo -e "sample_id\\tread\\treads_in\\treads_preumi\\tpct_kept_preumi"
    
    PCT_R1=\$(awk -v a="\${R1_IN}" -v b="\${R1_PRE}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
    printf "%s\\tR1\\t%s\\t%s\\t%.2f\\n" "\${SAMPLE_ID}" "\${R1_IN}" "\${R1_PRE}" "\${PCT_R1}"
    
    if [[ "\${MODE}" == "PE" ]]; then
      PCT_R2=\$(awk -v a="\${R2_IN}" -v b="\${R2_PRE}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
      printf "%s\\tR2\\t%s\\t%s\\t%.2f\\n" "\${SAMPLE_ID}" "\${R2_IN}" "\${R2_PRE}" "\${PCT_R2}"
    fi
  } > trim_stats.tsv

  echo "PREP | STATS | Trimming statistics saved: trim_stats.tsv"

  # UMI statistics (if UMI enabled)
  if [[ \${UMI_ENABLED} -eq 1 && \${UMI_LENGTH} -gt 0 ]]; then
    {
      echo -e "sample_id\\tumi_enabled\\tumi_source_read\\tumi_length\\tumi_location\\tumi_detect_mode\\tread\\treads_preumi\\treads_final\\tpct_kept_final"

      PCT_R1_UMI=\$(awk -v a="\${R1_PRE}" -v b="\${R1_FINAL}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
      printf "%s\\t1\\t%s\\t%s\\t%s\\t%s\\tR1\\t%s\\t%s\\t%.2f\\n" "\${SAMPLE_ID}" "\${UMI_READ}" "\${UMI_LENGTH}" "\${UMI_LOCATION}" "\${UMI_DETECT_MODE}" "\${R1_PRE}" "\${R1_FINAL}" "\${PCT_R1_UMI}"

      if [[ "\${MODE}" == "PE" ]]; then
        PCT_R2_UMI=\$(awk -v a="\${R2_PRE}" -v b="\${R2_FINAL}" 'BEGIN{print (a>0)?(b*100.0/a):0}')
        printf "%s\\t1\\t%s\\t%s\\t%s\\t%s\\tR2\\t%s\\t%s\\t%.2f\\n" "\${SAMPLE_ID}" "\${UMI_READ}" "\${UMI_LENGTH}" "\${UMI_LOCATION}" "\${UMI_DETECT_MODE}" "\${R2_PRE}" "\${R2_FINAL}" "\${PCT_R2_UMI}"
      fi
    } > umi_stats.tsv
    
    echo "PREP | STATS | UMI statistics saved: umi_stats.tsv"
  fi

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "PREP | README | Creating documentation..."

  cat > README_01_trimmed_fastq.txt <<'DOCEOF'
================================================================================
FASTQ PREPROCESSING — ${sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  This directory contains preprocessed FASTQ files ready for alignment.
  
  Processing Pipeline:
    Raw FASTQ → Cutadapt (adapters + barcodes + length) → UMI extraction → Final

FILES
────────────────────────────────────────────────────────────────────────────
  final_R1.fastq              — Cleaned R1 reads (ready for alignment)
  final_R2.fastq              — Cleaned R2 reads (or empty stub for SE)
  
  fastqc_raw/                 — Quality control reports on raw reads
  fastqc_final/               — Quality control reports on final reads
  
  trim_stats.tsv              — Read counts through trimming stages
  umi_stats.tsv               — UMI extraction statistics (if enabled)
  
  cutadapt_report.txt         — Detailed cutadapt trimming report
  umi_extract.log             — UMI extraction log (if enabled)
  preprocess_reads.log        — Complete processing log

PARAMETERS
────────────────────────────────────────────────────────────────────────────
  Sample ID:                  ${sample_id}
  Mode:                       ${(data_type ?: "SE").toString()}
  Threads:                    ${task.cpus}
  
  Adapter Trimming:           ${params.adapter_trimming?.enabled == true ? "enabled" : "disabled"}
  Barcode 1:                  ${params.barcode?.enabled == true ? "enabled (as configured: " + (params.barcode?.read ?: "R1") + " " + (params.barcode?.location ?: "5") + "', mode=" + (params.barcode?.detect_mode ?: "trust") + ")" : "disabled"}
  Barcode 2:                  ${params.barcode?.enabled2 == true ? "enabled (as configured: " + (params.barcode?.read2 ?: "R2") + " " + (params.barcode?.location2 ?: "5") + "', mode=" + (params.barcode?.detect_mode2 ?: "trust") + ")" : "disabled"}
  UMI Extraction:             ${params.umi?.enabled == true ? "enabled (as configured: " + (params.umi?.read ?: "R1") + " " + (params.umi?.location ?: "5") + "', mode=" + (params.umi?.detect_mode ?: "trust") + ")" : "disabled"}
  QC Tool:                    FastQC
  QC Enabled:                 ${(params.qc?.enabled == null ? true : params.qc?.enabled) ? "yes" : "no"}

  NOTE: "as configured" above reflects params.yaml, evaluated at pipeline
  submission time. If detect_mode is "auto" or "verify" for any slot, the
  VALUES ACTUALLY USED for trimming may have been overridden by the QC
  scan at run time -- see barcode_umi_profile.json and
  barcode_umi_detect.log (present in this folder when that happened) and
  the "PREP | DETECT |" lines in preprocess_reads.log for what was
  actually applied to these reads.

QUALITY METRICS
────────────────────────────────────────────────────────────────────────────
  See trim_stats.tsv for detailed read retention rates through each step.
  QC reports (HTML/JSON) are available in fastqc_raw/ and fastqc_final/.

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  Use final_R1.fastq (and final_R2.fastq for PE) for alignment.
  These files have been:
    ✓ Adapter trimmed (if enabled)
    ✓ Barcode removed (if enabled)
    ✓ UMI extracted and appended to read names (if enabled)
    ✓ Length filtered (minimum insert length enforced)
    ✓ Quality checked

NOTES
────────────────────────────────────────────────────────────────────────────
  • Single-pass processing: Faster than traditional multi-step approaches
  • UMI information (if extracted) is embedded in FASTQ read names
  • Empty final_R2.fastq for SE samples maintains consistent file structure

================================================================================
DOCEOF

  echo "PREP | README | Documentation created"

  ###########################################################################
  # 10) SUMMARY
  ###########################################################################

  # Get final file sizes
  FINAL_R1_SIZE=\$(stat -c%s final_R1.fastq 2>/dev/null || stat -f%z final_R1.fastq 2>/dev/null || echo "unknown")
  FINAL_R2_SIZE=\$(stat -c%s final_R2.fastq 2>/dev/null || stat -f%z final_R2.fastq 2>/dev/null || echo "unknown")

  echo "────────────────────────────────────────────────────────────────────────"
  echo "PREP | SUMMARY | Processing complete for \${SAMPLE_ID}"
  echo "PREP | SUMMARY | Mode: \${MODE}"
  echo "PREP | SUMMARY | Input reads (R1): \${R1_IN}"
  if [[ "\${MODE}" == "PE" ]]; then
    echo "PREP | SUMMARY | Input reads (R2): \${R2_IN}"
  fi
  echo "PREP | SUMMARY | Final reads (R1): \${R1_FINAL} (\$(awk -v a="\${R1_IN}" -v b="\${R1_FINAL}" 'BEGIN{printf "%.1f", (a>0)?(b*100.0/a):0}')%)"
  if [[ "\${MODE}" == "PE" ]]; then
    echo "PREP | SUMMARY | Final reads (R2): \${R2_FINAL} (\$(awk -v a="\${R2_IN}" -v b="\${R2_FINAL}" 'BEGIN{printf "%.1f", (a>0)?(b*100.0/a):0}')%)"
  fi
  echo "PREP | SUMMARY | Final R1 size: \${FINAL_R1_SIZE} bytes"
  echo "PREP | SUMMARY | Final R2 size: \${FINAL_R2_SIZE} bytes"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "PREP | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}