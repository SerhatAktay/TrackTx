# Changelog

All notable changes to TrackTx are documented in this file.

## [Unreleased]

## [1.4.0] - 2026-09-17

54 commits since 1.3.0: correctness fixes across the compute path (divergent-transcription calling, pol metrics, normalization), a new barcode/UMI auto-detection mode, first nf-test coverage, and a full README/CHANGELOG restructure.

> Note: this repo also has an older tag `v1.5` (2025-12-19), from before the Nextflow 26 rewrite that reset numbering back through v1.2.0 → v1.3.0. It predates and is unrelated to this 1.4.0 line; ignore it when comparing tags by semver order.

### Added
- New barcode/UMI auto-detection mode (`bin/detect_barcode_umi.py`, wired through `modules/03_preprocess_and_quality_filter_reads.nf` and the config generator) so a sample's barcode/UMI layout no longer has to be specified by hand.
- `modules/00_capture_tool_versions.nf`: dumps the shared conda/micromamba environment's resolved package versions plus pipeline/Nextflow version once per run to `pipeline_info/versions.yml`, giving a single citable provenance file instead of grepping `.nextflow.log`.
- `nextflow_schema.json` + nf-schema plugin: validates the ~25 flags a user actually types (required fields, enums) on launch; existing hand-written conditional validation in `main.nf` remains the source of truth for logic the schema can't express.
- Empirical-null FDR check for the divergent-transcription caller: shuffles each strand's signal within-chromosome and reports a measured false-positive rate alongside the GMM-posterior score in the QC report (see Fixed/Changed below for what this uncovered).
- First nf-test coverage (modules 01, 04, 08, 09, 10, 11, 12) plus `bin/*.py` self-tests, wired into `.github/workflows/ci.yml` alongside the existing `-preview` structural-only job -- `-preview` never executes a process body, so a broken command or bad flag previously passed CI silently.
- Eight real-dataset example `params/*.yaml` + `samplesheets/*.csv` pairs (human, mouse, dog, drosophila, arabidopsis) as ready-to-adapt starting configs.
- README split into `docs/MODULES.md`, `docs/INSTALLATION.md`, `docs/TROUBLESHOOTING.md`, trimmed and de-duplicated in the process.

### Fixed
- `bin/detect_divergent_transcription.py`: the GMM's positive component was picked by mean signal (`log_total`), not by strand balance -- since "divergent" is defined by bidirectional balance, not read depth, a high-signal one-sided region (readthrough, alignment artifact) could outrank a genuinely balanced but lower-signal one. Component selection now runs on `balance_bayesian` first, with `log_total` demoted to a logged cross-check; changes ~1/3 of calls on real data (2,015 dropped, 2,296 gained on a K562 test sample) and lifts mean balance among passing calls from 0.34 to 0.88. Also fixed `merge_overlapping_regions()` taking `max()` instead of summing the `total` signal column across merged peak-pairs, undercounting merged regions' reported signal.
- Investigated the caller's empirical FDR (measured ~0.51-0.53 against a nominal 0.05 cutoff) as a possible calibration bug; four independent approaches (extra GMM feature, local-background threshold, raw stringency, target-decoy cutoff calibration) failed to close the gap without collapsing the catalog to near-zero calls. Concluded this is a real ceiling of peak-pair-aggregate scoring against pervasively-transcribed PRO-seq signal, not a solvable tuning problem -- every place that implied a calibrated FDR (module header, embedded docs, `--fdr` help text, QC report wording) now says plainly this is score-based ranking with a measured ~0.5 real error rate, pointing to the QC report's own measurement each run.
- `modules/05b_check_and_merge_replicates.nf`: the replicate-concordance check computed correlation over per-chromosome total read counts (`samtools idxstats`), dominated by chromosome length and with little power to catch real discordance. Replaced with genome-wide binned read counts (`bedtools makewindows` + `bedtools multicov`, default 10 kb bins, `params.replicates.concordance_bin_size`), with counts log2(n+1)-transformed before Pearson (raw Pearson over-penalized ordinary noise at low/moderate counts -- pre-fix output read 0.58-0.85 against 0.93-0.98 in source papers for the same data; log2 transform reproduced published values within ~0.01-0.02) and zero-count bins no longer excluded (the previous `$c > 0` filter was silently discarding exactly the bins that most reveal discordance).
- `bin/calculate_pol_metrics.py`: TSS/body reads now counted by 3'-end (Pol II active-site) position instead of full-read overlap, matching the TSN convention used elsewhere and removing a read-length-dependent bias in the pausing index; `required_body_len` now folds in the per-organism `auto_body_offset` calibration as a floor instead of staying mammalian-scale for every organism; gene/neighbor windows are clipped at the overlap midpoint with the nearest same-strand neighbor so compact genomes (bacterial operons, dense plant/insect genomes) don't bleed adjacent-gene signal into the pausing index.
- `modules/11` & `modules/13`: `samtools` only honors the last `-F` flag when given twice; combining the exclude flags into one value fixes the unmapped-read exclusion silently being dropped whenever dedup was enabled (the default).
- `download_genome_annotations.nf`: the genome-cache fast-path cache key was missing `annotation_source`/`annotation_exclude_biotypes`/`annotation_chr_naming`, so a rerun with a different annotation config could silently reuse a prior config's cached GTF/gene-catalog files.
- Shared ERR trap called its error helper with only 3 args, always exiting 1 regardless of the failing command's real exit status -- silently masked OOM kills (137/143) as a generic failure, defeating every process's retry-with-more-memory `errorStrategy`.
- `run_batch.sh` (now removed, see Removed) deleted `work/.nextflow` state after every job including failed ones, defeating `-resume` for exactly the case that mattered.
- Module 09 `SCORE_MIN`/`SCORE_MAX` crashed under `pipefail` on a `sort | head` SIGPIPE; module 03's barcode/UMI verify-mode error path had a bash syntax error; module 1 had a bug when using a custom/"other" reference genome; module 1's README heredocs across modules 02/05/06/07/08/09/10/11/12 had quoted delimiters so runtime values (counts, dates, stats) printed as literal `${VAR}` text instead of substituting.
- `TrackTx_config_generator.html`: `buildCSV()` had no CSV escaping (a comma in a free-text field silently shifted every later column); YAML string-quoting escaped `"` but not `\` (a Windows-style path corrupted or hard-failed YAML parsing); the custom-genome GTF requirement check used a bare `alert()` with no `return`, so it never actually blocked.
- `nfmon.py`: `tail_trace()` advanced its read pointer past an incomplete trailing line before checking it parsed, permanently dropping that line's status update if Nextflow was mid-write when polled (a task could show stuck "RUNNING" for the rest of the session); its Run/Session-name regexes never matched Nextflow 26's actual log wording, so both fields stayed at "?" on every run; added a SIGTERM handler (only SIGINT was caught, so a plain `kill` left the terminal in raw curses mode); narrowed the Rich-UI-unavailable fallback to `ImportError` only, so a real bug in the Rich path no longer silently launders into curses mode. Several earlier passes this cycle fixed additional nfmon UI/display bugs and retuned module 05b's CPU allowance.
- `run_pipeline.sh`: disk-space check always checked `.` even in `--external-drive` mode, where the real work/cache destination is `${HOME}/tmp`; added a profile-consistency check on `-resume` (a silent profile mismatch risks a full recompute since Nextflow's cache hash includes container/conda directives).
- `storeDir` migration (see Changed) exposed colliding `3p/**`/`5p/**` glob outputs in `modules/08_normalize_coverage_tracks.nf` that matched files already claimed by named outputs earlier in the block; both emits were unused by any consumer and removed.
- Fixed sci-notation CPM values silently breaking `bedGraphToBigWig` into empty tracks, malformed sample files silently merging as NaN rows, one adapter-dimer read collapsing barcode/UMI detection, PE fragments double-counted from missing mate filtering, and genome-size auto-scaling misclassifying real genomes from sparse bedgraph coordinates alone.
- Swept all 7 `datetime.datetime.utcnow()` call sites (deprecated) across 5 scripts to `datetime.now(timezone.utc)`, with the isoformat-plus-literal-`"Z"` sites additionally fixed to avoid a malformed doubled `"+00:00Z"` suffix; output format unchanged.
- `_common.py`'s new `AtomicFileWriter` (see Changed) initially left every atomically-written output owner-only-readable (mode 0600, inherited from `tempfile.mkstemp` via `os.replace()`); fixed to chmod before writing.
- si-normalization silently not applying as intended, and a gene-list bug in the GTF-to-catalog script caused by GTF file structure, both traced and fixed; `control_label`'s default (`'CTRL'`) is now documented in `nextflow.config` as a placeholder that must be overridden per dataset, not a usable fallback.

### Changed
- Unified Python logging, `__main__` error handling, and output writes across all 8 `bin/*.py` scripts via new `bin/_common.py` (previously 3 scripts duplicated identical logging helpers, one used a different logger, three had none at all; 4 different `__main__` error-handling patterns existed with no majority). All primary-output writes now go through `AtomicFileWriter` (temp-file-in-same-dir + `os.replace()`), matching the mktemp+mv pattern already used on the shell side, so a killed process can no longer leave a truncated file at the real output path.
- Consolidated per-file version-history narration that had accumulated directly in source comments (`bin/functional_regions.py`'s v9.0/v9.1 rewrite notes, `bin/calculate_pol_metrics.py`'s `auto_body_offset` tuning history) into this changelog; source comments now state current behavior plus a short rationale and point here for history.
- Unified header style across `bin/*.py` (previously 3 different conventions: ASCII banner, Unicode box-drawing, bare docstring) and made `modules/11b_collect_pol_metrics_per_replicate.nf`'s header match its 15 sibling modules; removed `detect_divergent_transcription.py`'s per-file `Author`/`Date` byline (authorship belongs in repo metadata, not scattered per-file).
- `publishDir` switched to `storeDir` across modules 03/05/05b/06/08, then which modules keep persistent storage under it was re-scoped in a follow-up pass.
- Auto-scaled the divergent-transcription caller's `nt_window`/`bin_gap`/`merge_gap`/`bg_window` from apparent genome size (same pattern as `auto_body_offset`) instead of fixed human-sized defaults, so a compact genome's window doesn't span multiple adjacent genes; fixed a latent mismatch this surfaced where `main.nf`'s default disagreed with the script's own argparse default.
- Deduplicated the custom-genome-id sanitize regex (previously reimplemented independently in three places) into `lib/GenomeId.groovy`, shared by `main.nf` and `download_genome_annotations.nf`.
- Retuned alignment/divergent-detection CPU tiers and module 05b's CPU allowance; threaded `THREADS`/`-@` consistently through modules 05/06/07/08/10/11/13's `samtools`/`sort` calls (previously inconsistent, some effectively single-threaded); parallelized module 08's normalize step and UMI dedup; vectorized hot loops in `detect_divergent_transcription.py` and `functional_regions.py`; removed dead config flags (`qc.depth_max_cov`, `preprocess_reads_lenient_cache`) and dead script args (`count_mode`, `div_fallback_*` in `functional_regions.py`) uncovered while retuning.
- Module 02's cache check now requires non-empty content plus a FASTQ header sanity check for uncompressed files, not just file existence, so a truncated file left by an interrupted prior run is no longer trusted as "already downloaded."
- README restructured and trimmed for clarity (split into `docs/`, see Added); documented why TrackTx ships one container/conda env for every process instead of nf-core's per-tool-container convention.

### Removed
- `scripts/run_batch.sh` and its example job template, per explicit request -- superseded by other tooling; see Fixed for the `-resume`-defeating bugs this removal makes moot.
- Code-review working notes untracked from the repo; `params/`, `samplesheets/`, and other run-local artifacts untracked via `.gitignore`.

### Notes for future work (flagged, not changed)
- `bin/functional_regions.py`'s sequential read-assignment step removes every overlapping read from further consideration on an *unstranded* basis at each stage (Promoter, Divergent, CPS, Gene body, Termination), even though assignment itself is strand-specific. This means an antisense read overlapping a region's window is dropped before it can be assigned to that category, reach a later category, or land in "Non-localized polymerase" -- a real but currently unreported source of read loss. This matches the original TrackTx.sh logic intentionally (see `sequential_read_assignment()`'s docstring for the mechanism); flagging here since it affects how much of total input reads are accounted for in `functional_regions_summary.tsv`.
- The divergent-transcription caller's empirical FDR ceiling (~0.5, see Fixed) is reported but not gated on; a shape-based redesign closer to how dREG/PINTS work internally is out of scope for now.

## [1.3.0] - 2026-07-01

This release promotes the `dev` branch (Nextflow 26 compatibility line) to `main`.
It does **not** include a separate line of work that had been developing on the old
`main` in parallel (GRO-seq support, native TAIR10 handling, TSS-based replicate
merging, an nfmon TUI redesign, and additional align/divergent-detector performance
tuning). That history is preserved under the tag `main-archive-pre-v1.3.0` for future
porting into `dev`/`main` on a case-by-case basis.

### Added
- Multimapper-aware `allMap` coverage tracks alongside the existing `main` tracks, with `NH`-tag-based uniqueness (`align.multimap_k`) replacing MAPQ-only filtering throughout QC, pol metrics, and reports.
- New cohort QC & visualization module: MultiQC aggregation, deepTools PCA/correlation heatmap on CPM-normalized tracks, an IGV session file, and a run-on efficiency table.
- `sra_cache_dir` param and additional results-size controls (toggle bedGraph/allMap/5′-track publishing) to shrink output footprint.
- Memory-derived concurrency cap for coverage-track generation to prevent OOM when multiple `genomecov` jobs run at once.

### Changed
- Full Nextflow 26 compatibility rewrite across `main.nf` and all modules (strict config-parser handling, inline resource-allocation helpers since `lib/*.groovy` classes are no longer visible to config closures).
- SRA/FASTQ downloads now build their channel from `storeDir` cache on disk, fixing a deadlock that occurred once all samples were already cached.
- Module 01 genome/annotation downloads: ranged-GET UCSC probe (was a broken HEAD request) with a UCSC-rsync → NCBI-RefSeq fallback, including chromosome-name normalization.
- siCPM normalization channel wiring fixed, unlocking `output.bedgraph: false` as a supported combination.
- Pol-aggregate step now uses a correct t=0 baseline; differential contrasts are force-disabled by design when replicates are merged (superseded by the log2FC approach).

### Fixed
- CRLF line endings from NCBI assembly reports leaving stray `\r` in chromosome names, which crashed module 11 (pol metrics) and caused `bedtools` OOM.
- Divergent-transcription detector locale bug and heatmap/report timepoint columns sorting lexicographically instead of numerically (e.g. `0, 10, 160, 20, 40, 60`).
- Various module-specific OOM and bug fixes (modules 01, 05, 06, 08, 11, 13, 14) surfaced during Nextflow 26 migration and full-cohort validation runs.

### Performance
- Vectorized the divergent-transcription detector (prefix-sum `FastBedGraph`, dropped per-row `iterrows`) — cut per-sample runtime from ~2 hours to a fraction of that, output validated byte-identical to the previous implementation.

### Docker
- Default pipeline image tag bumped to `ghcr.io/serhataktay/tracktx:1.3.0` (`manifest.version` in `nextflow.config`). **This image must be built and pushed to GHCR before running `-profile docker`/`singularity`/`podman` against this release** — see the release instructions for the build/push commands.
