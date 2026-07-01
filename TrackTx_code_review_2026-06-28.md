# TrackTx — full code review (data-flow audit)

> **STATUS — fixes applied 2026-06-28 (dev repo `SciLifeLab/tracktx` only; SanDisk run copy NOT touched).**
> Resolved: **A1** (run-on now parses the genes.tsv schema), **A2** (cohort lists projected from one
> sorted tuple so conditions stay aligned), **A3** (siCPM control picks lowest replicate incl. merged
> rep-0), **A4** (report reads an existing QC-JSON key), **C2** (PE Pol-II gene metrics now use the same
> mate-filtered BAM as the tracks), **D3** (5′ auto-detect uses `-s`), **F3** (download processes retry),
> **F4** (dead commented validations removed), **G** (pausing-index validator 7→8), **E** (dead files
> deleted: 10a, enhancer_gene_score.py, normalize_tracks.py, .bak/.deprecated, __pycache__; dead
> `_customAnnotationFile`/NO_GTF block removed), **H** (divergent + module-06 doc drift corrected).
> Decisions you made (round 1): **B1** → median-gated body offset (≥10 kb median ⇒ keep user 2000 bp; only
> compact genomes shrink); **B2** → kept GMM caller, the top-fraction fallback is now **opt-in**
> (`advanced.divergent_fallback_top_frac`, default 0 ⇒ a noisy sample legitimately returns 0 sites) and
> the `--fdr`/README wording clarifies it's an approximate posterior cutoff; **B4** → kept host-first spike
> alignment (documented).
>
> **Round 2 (also done):** **B3** → Pol-II metrics now read the gtf_to_catalog `genes.tsv` (new
> `parse_catalog_file` + `--genes`; module 11 + main.nf wired to `genes_ch`) so the pausing-index TSS/TES
> are identical to functional-region calling. **G** → `functional_regions.py` matches active genes by unique
> `gene_id` instead of symbol (no more paralog collisions). **F1** → per-sample report track links are now
> derived from the REAL produced files on channels (module 08 `report_bw` emit + module 06 allMap raw
> bedGraphs), gated by publish toggles, with the canonical published path as link text — this also fixed a
> latent wrong `cpm/3p/` publish path that never resolved; `resolveOutputPath` removed. **B2 follow-on:**
> module 10 + the detector now treat an EMPTY divergent BED as a legitimate 0-site result (warn + continue)
> instead of a fatal error, so "return zero sites" works end-to-end.
>
> **Round 3 (final pass — remaining items):** **C1** → the PE signal-mate is now a parameter
> (`align.pe_signal_mate`, default `'read2'` = flag 128 = current behavior; `'read1'` = flag 64 for other
> chemistries), replacing the hard-coded `-f 128`. **C3** → optional allMap UMI-dedup
> (`umi.dedup_allmap`, default false = current behavior) so allMap tracks can match the deduped main tracks;
> off by default with a note that umi_tools on -k multimappers is approximate, and it falls back to the raw
> allMap BAM on failure. **D4** → module 07 header now states plainly that it computes whole-library totals,
> not per-gene counts (the process can't be renamed without breaking the config `withName` selectors).
> **C4** → re-examined and found to be a NON-issue (genome and spike counts are both per-read-record, so the
> PE spike-fraction is dimensionally consistent); annotated so it isn't "fixed" wrongly. **F2** (hand-rolled
> SRR resume) and **F5** (PE un-conc spike feed) left as-is — both work and rewriting them is higher risk
> than value; D1/D2 are intentional cost trade-offs.
>
> Verification (all rounds): `py_compile` on every `bin/` script + every embedded heredoc, `bash -n` on the
> rendered script body of all edited modules, brace balance on all edited `.nf` + config.
>
> **LIVE NF26 VALIDATION (2026-06-29):** stood up Temurin JDK17 + Nextflow 26.04.3 in-sandbox and ran the
> real two-layer check. This surfaced and fixed **two parse-blocking errors that would have crashed at
> launch**: (1) a pre-existing Java-style `(double) gmap` cast in main.nf STEP 6c (NF26 strict parser
> rejects it → rewritten `(gmap as double)`); (2) my F1 `has(x)` named-local-closure call (same strict
> gotcha → typed `Closure<Boolean> has` + `has.call(x)`). After the fixes: **`nextflow lint` = 0 errors**
> (only advisory `projectDir` warnings); **`nextflow run -preview` (full DAG) exits 0** under the default
> config, `--replicates.merge true`, and the toggled combo (`umi.dedup_allmap`, `align.pe_signal_mate=read1`,
> `output.bedgraph=false`, `norm.emit_allmap=false`); **`nextflow config` resolves** on
> docker / conda,performance / slurm,singularity. NOT executed: real tasks (need Docker + the ~3 GB hs1
> genome, which won't fit the sandbox), so task-internal bash/python runtime behavior is still only validated
> by `bash -n`/`py_compile` — confirm with a real `-profile docker` test_PE run on your Mac.

> **Round 4 — real test_PE run analysed (2026-06-29).** A live `-profile docker` test_PE run (1 PE K562
> control sample, T2T-hs1 + dm6 spike) completed all 12 stages with no `tracktx_error`. Confirmed working
> from real outputs: **A1** run-on now uses 13,926 gene bodies (was 0); **A3** siCPM factor non-zero and
> equals CPM for the control (correct by design); **B2** 5,848 divergent regions; **B3** pol-metrics ran on
> the catalog (57,514 genes, `Gene source: catalog hs1.genes.tsv`); **G** functional regions populated
> (3,472 active genes). The run ALSO surfaced two more bugs, now fixed: **B1** — the median gate
> mis-classified human/T2T as "compact" (all-gene catalog median ~5.9 kb < the 10 kb threshold) and shrank
> the offset to 200 bp; threshold lowered to 4 kb so human/mouse keep the user 2000 bp. **F1** — links were
> computed correctly but module 14 re-`stat`'d them against its task sandbox (publish dir not visible there)
> and dropped all 6 ("0/6 tracks"); module 14 now trusts the upstream-validated non-empty link and main.nf
> emits absolute paths. Two non-code observations to revisit on full data: the run-on 5′/3′ metric clusters
> at ~1.0 (the two ends of the same short PRO-seq read both fall in the body, so the ratio is near-tautological
> — metric DESIGN worth rethinking), and the test's pausing-index magnitudes (median ~174, 74% truncated) are
> dominated by 10%-subset sparsity, not validatable here.

> **Round 5 — the two flags fixed (2026-06-29).** **Run-on metric** redefined from the near-tautological
> 5′-end-track ÷ 3′-end-track ratio to a POSITIONAL distal(3′)/proximal(5′) Pol II (3′-end) signal ratio over
> long gene bodies — the real 5′→3′ falloff measure, and it now uses only the 3′ track so it also works for
> single-end (the old one was skipped for SE). Validated on the real test 3′ bedGraphs: **median = 0.52
> ("good") over 7,776 genes** vs the old 1.0000. Landing-page wording updated. **Pausing-index summary**: the
> per-sample report median PI now uses a body-count coverage floor (`pol.pi_min_body_count`, default 10;
> falls back to all valid genes if <20 clear it), so the headline isn't inflated by 1-read genes on shallow
> data; the per-gene `pausing_index.tsv` is unchanged and the JSON records `pausing_genes_used`. Re-verified:
> lint 0 errors, `-preview` exit 0 (default + `emit_5p=false`), all heredocs `py_compile`.

Date: 2026-06-28. Scope: `main.nf`, `nextflow.config`, all 16 modules, all `bin/` scripts.
Method: traced one sample (and a cohort) end-to-end through every process and python/awk script,
looking for places where the data is mishandled, where behaviour silently deviates from intent,
where work is done twice, and where logic only holds for one specific kind of input.

The pipeline is, on the whole, coherent and a real PRO-seq pipeline. But a lot of recent
bug-fixing has been *local* — each module was patched in isolation — and the **seams between
modules** are where most of the remaining problems live. Several are silent: the run still
reports success and produces plausible-looking files.

Issues are ordered by how much they distort results, not by how hard they are to fix.

---

## A. Correctness bugs that silently produce wrong / empty results

### A1. Run-on efficiency is fed the wrong file → metric never works *(high)*
`main.nf` passes `genes_ch` (= `…​.genes.tsv`, the 9-column catalog
`gene_id, gene_name, chr, strand, start, end, tss, tes, biotype`) into
`cohort_qc_and_viz` as `genes_bed`. The run-on python (`16_…nf`, lines ~474-501) parses it as
**BED6** and guards each row with:

```python
if not (p[1].lstrip('-').isdigit() and p[2].lstrip('-').isdigit()): continue
chrom, start, end, name, _, strand = p[0], int(p[1]), int(p[2]), p[3], p[4], p[5]
```

In `genes.tsv`, `p[1]` is `gene_name` (not a number) and `p[2]` is `chr` (not a number), so
**every gene row is skipped**. The result is always 0 gene bodies → `runon_efficiency.tsv` is
always "insufficient_data" (or, for SE, the not-applicable path). The landing-page Run-on table
is therefore always empty/meaningless. Fix: pass `tss_bed`/`tes_bed` (real BED6) or a derived BED6, or
teach the parser the genes.tsv column order. This is the single most clearly-broken feature.

### A2. Cohort lists are sorted independently → condition labels can mis-map to samples *(high)*
`main.nf` STEP 16 builds each cohort list with its own `.toSortedList()`:

```
cohort_bw_pos3   …→ bwp.toString() }.toSortedList()
cohort_sample_ids…→ sid           }.toSortedList()
cohort_conditions…→ c ?: 'unknown'}.toSortedList()
```

The BigWig/bedGraph/sample_id lists happen to co-sort because every BigWig path embeds the sid,
so element *i* lines up across those. **`cohort_conditions` does not** — it is sorted by the
condition string itself. So in `cohort_qc_and_viz`, `CONDITIONS[i]` is not guaranteed to be the
condition of `SAMPLE_IDS[i]`. This corrupts the IGV per-condition colour grouping and any
condition-based interpretation of the deepTools panel. Fix: carry one joined tuple and emit the
columns from a *single* `toSortedList` (sort once by sid, then map out each field), instead of
sorting each field separately.

### A3. siCPM control selection breaks when replicates are merged *(high, if spike-in used)*
`08_normalize`'s control picker looks for `condition == control_label AND replicate ∈ {1,r1,R1}`.
But merged samples are emitted with **replicate = 0** (`05b`/`main.nf`), so the control's merged
track never matches, and it silently falls through to "first sample with spike_reads>0" — i.e.
the first row alphabetically, which is usually the *wrong* sample. With `replicates.merge=true`
+ spike-in, siCPM is therefore scaled to an arbitrary reference. Fix: match `replicate in {0,1,…}`
or select the control by condition only (then by smallest replicate).

### A4. Per-sample report parses a QC-JSON key that doesn't exist *(low, cosmetic)*
`14_…nf` reads `jq -r '.map_rate_percent'` but `13_…nf` writes `primary_retained_percent` /
`mapq_pass_percent` (there is no `map_rate_percent`). So `MAP_RATE` in the per-sample README is
always 0/NA. The true overall alignment rate is also never threaded into the per-sample HTML
(by design it lives only in `02_alignments/alignment_rates_summary.tsv`). Net: the per-sample
report cannot show a real mapping rate.

---

## B. Scientific / methodological deviations from intent

### B1. Pol-II "auto body offset" silently shrinks the body window for *every* organism *(medium)*
`calculate_pol_metrics.py::auto_body_offset` returns `min(user_offset_min, max(200, P25*0.20))`.
The intent (per docstring) was to *reduce* the offset only for compact genomes. But P25 is taken
over **all transcript/exon-feature lengths**, which is small even in human/mouse (lots of short
transcripts), so `P25*0.2` is routinely < the 2000 bp default. Because it takes the `min`, the
human body window starts at e.g. TSS+1.0–1.4 kb instead of TSS+2 kb. That pulls the
promoter-proximal pause peak partly into the "body" and **deflates the pausing index
systematically**. If the goal was "never *larger* than user, but only shrink for small genomes",
the calibration should compare against the *gene-length median*, not P25 of all features, and
probably gate on organism scale.

### B2. Divergent detector is a bespoke GMM heuristic, not a calibrated/【published】caller *(medium)*
`detect_divergent_transcription.py` does threshold→peak→pair→2-component GMM→"FDR". Several points
worth a decision:
- The "FDR" is a posterior-based approximation, **not** a real BH FDR (the code itself says so).
  The `divergent_fdr` knob therefore doesn't mean what users will assume.
- When nothing passes, it falls back to "top 10% or 100 regions" — guaranteeing output even from
  noise.
- With a single candidate pair it emits it at confidence 1.0.
- One-positive-to-many-negative pairing inflates raw pair counts, then `merge-gap` collapses
  them; the final count is sensitive to `merge_gap` and `nt_window` in non-obvious ways.
This is fine as an in-house method but it is not equivalent to the original TrackTx fixed-window
divergent definition, and the QC/README text still describes the *old* 95th-percentile / 10×
parameters (see H1). Decide whether this is the intended scientific method and document the real
behaviour.

### B3. Two different gene models feed two different steps *(medium)*
- Functional-region calling (`10`) and run-on use the **`gtf_to_catalog` genes.tsv** (one
  TSS/TES per gene, gene-feature-derived).
- Pol-II metrics (`11`/`calculate_pol_metrics.py`) **re-parse the raw GTF** and pick the *longest
  transcript per gene*.
So the TSS used for "promoter signal" (module 10) and the TSS used for the pausing index
(module 11) can differ for the same gene. For internal consistency both should derive from the
same catalog. At minimum, document that promoter-region geometry and pausing-index geometry use
different TSS definitions.

### B4. Spike-in is aligned only to genome-*unmapped* reads *(medium — a real design choice)*
`05_align` sends only `unaligned.fastq` (reads that failed the primary genome) to the spike
genome. This is the "competitive/host-first" convention and is defensible, but it means any read
that maps to *both* host and spike is assigned to host, so `spike_mapped_reads` (the siCPM
denominator) is a lower bound and depends on host mapping stringency. Standard spike-in workflows
often align all reads to a combined index. Worth a conscious decision since siCPM scaling rides
entirely on this count.

### B5. Functional-region signal uses CPM tracks; region *definition* uses raw tracks *(low)*
Module 10 assigns **raw** 3′ signal to regions (good, matches old script). Module 11 density then
re-quantifies the same regions from **CPM/siCPM** tracks. Two different signal bases for "signal
per region" land in different output files (`functional_regions_summary.tsv` vs `pol_density.tsv`)
with no flag telling the reader they aren't comparable. Not wrong, but a trap for downstream users.

---

## C. Paired-end assumptions (only correct for one library design)

### C1. The whole PE path hard-codes one specific PRO-seq mate layout *(medium)*
`05_align` does `bowtie2 --ff -1 R2 -2 RC(R1)` and then everything downstream assumes
**flag-128 (second-in-pair) = RC(R1) = the Pol II 3′ position**:
- `06_tracks` keeps only `-f 128` reads for coverage.
- The comments encode this as fact.
This is correct *only* for that exact chemistry/adapter orientation. There is no parameter to
describe a different PE layout; a user with standard PE PRO-seq of the opposite orientation gets
silently wrong tracks. Consider deriving mate-of-interest from `library_type`/an explicit param
rather than baking `-f 128` in.

### C2. Pol-II gene metrics count *both* mates in PE *(medium)*
`11` feeds `bam_for_downstream.bam` (the deduped main BAM, **both mates**, not the `-f 128`
Read2-only BAM that tracks use). `calculate_pol_metrics.py` then counts every primary read whose
strand matches the gene. In PE this counts the R2/"wrong-end" mate too, so TSS and body counts are
inflated/contaminated relative to the coverage tracks. The pausing index partly cancels (ratio),
but `tss_cpm`/`body_cpm` and the per-replicate hand-off table do not. Tracks and gene-metrics
should use the *same* mate-filtered BAM in PE.

### C3. allMap BAM is never UMI-deduplicated *(low)*
In `06`, UMI dedup is applied to the main BAM only; allMap coverage (and allMap CPM tracks, and
the allMap tracks shown in IGV/divergent docs) come from the **non-deduplicated** allMap BAM. So
"main" and "allMap" tracks are deduped inconsistently when `umi.enabled=true`.

### C4. spike_fraction mixes PE mate-counted genome reads with SE-counted spike reads *(low)*
`main.nf` STEP 6c computes `spike_fraction = spike_mapped/genome_mapped*100`. In PE,
`genome_mapped` counts mate records (≈2× pairs) while spike alignment is SE (1× per read), so the
reported spike fraction is ~halved in PE. Cosmetic but misleading for the "is there really a
spike-in?" check.

---

## D. Redundant or doubled work

### D1. Pol-II gene metrics are computed twice when merging *(medium, cost)*
With `replicates.merge=true`, `calculate_polymerase_occupancy_metrics` runs on the merged BAM
**and** `pol_metrics_per_replicate` runs the *same* heavy calculation on every individual
replicate BAM (11b path). That's intentional (differential hand-off needs per-rep), but the
per-replicate run also recomputes the density/functional inputs it then throws away (passed empty
placeholders). The merged gene-metric pass is effectively redundant with the per-rep passes for
anything except the cohort table. Worth confirming both are actually needed.

### D2. bedGraph → BigWig and re-sorting happen in several places *(low, cost on USB)*
`06` sorts during `genomecov | sort`; `08` can re-sort (`force_sort_bedgraph`); `make_bigwig`
helpers exist in both `06` and `08` with subtly different memory logic. Multiple
write→read→rewrite passes of multi-GB bedGraphs is the dominant cost on the SanDisk run copy.
The 5′ tracks are *always* generated in `06` and *always* normalized in `08` (see D3) even when
nobody consumes them.

### D3. 5′ tracks are effectively always-on even when "auto" should skip them *(low)*
`08`'s `EMIT_5P` auto-detect is `if [[ -e "$POS5" || -e "$NEG5" ]]`. `main.nf` always supplies a
sentinel `EMPTY_5P_POS.bedgraph` that *exists*, so auto mode always enables 5′ normalization,
which then processes the (often empty) inputs and writes empty outputs. The intended "skip when no
5′ data" never triggers. Either test for non-empty (`-s`) or drive it from `library_type`.

### D4. Counts are collected by idxstats, never used as a gene-level quantification *(low)*
`07_quantify_reads_per_gene` is named "per gene" but only sums `idxstats` (whole-BAM main/allmap/
spike totals). That's all `08` needs (library sizes), but the name and the README promise gene
counts that don't exist. Rename or document.

---

## E. Dead / unwired code (confuses maintenance, invites accidental use)

- `modules/10a_score_enhancer_vs_gene.nf` + `bin/enhancer_gene_score.py`: **not imported** in
  `main.nf`, never run. (10a even re-declares `nextflow.enable.dsl = 2` mid-module.)
- `bin/normalize_tracks.py`: **no references anywhere** — normalization is done inline in awk in
  module 08. 370 lines of dead Python that looks authoritative.
- `bin/detect_divergent_tx_old.py.deprecated`, `bin/*.bak.*`, `bin/combine_reports.py.aibak`:
  stale copies sitting next to the live scripts.
- `main.nf` `_customAnnotationFile` (lines 226-231) is computed and never used;
  `download_genome_annotations()` is called with **no arguments**, so the workflow-level
  custom-GTF resolution is dead (the module still honours `params.gtf_path` internally, but only
  when `reference_genome == 'other'`).
- `__pycache__/` and a `.fuse_hidden…` file are checked into the working tree.

---

## F. Robustness / caching traps

### F1. Report inputs are read back from the *publish* directory, not channels *(medium)*
`main.nf` STEP 15 builds track links with `resolveOutputPath("${params.output_dir}/…")` — i.e. it
reaches into the published results layout and tests file existence at channel-map time. This
couples the report to (a) publish having completed, (b) `publish_mode`, and (c) every
space-saving toggle (`output.bedgraph=false`, `emit_allmap=false`, …). When those are off the
links silently become empty. The robust pattern (used elsewhere now) is to pass the files through
Nextflow channels, not to re-derive publish paths.

### F2. SRR resume logic is a hand-rolled disk-cache state machine *(medium)*
STEP 3's "trimmed exists? cached raw exists? else download" branching (plus the `storeDir`
deadlock history in memory) reimplements caching that Nextflow's own `-resume`/`storeDir` is meant
to handle. It works, but it's fragile: the `sraCacheDir` string in `main.nf` must stay byte-for-
byte in sync with the `storeDir` closure in `02`, and `by_trimmed`/`by_cache` evaluate
`file().exists()` at runtime against the results drive. Any change to `publish_sra_fastq` /
`sra_cache_dir` defaults can re-trigger downloads or deadlocks.

### F3. `errorStrategy` only retries OOM (137/143); everything else terminates the run *(low)*
Transient network failures in download modules (`01`,`02`,`04`) are not retried at the Nextflow
level (only in-script `curl --retry`). A single flaky UCSC/NCBI/ENA hiccup outside curl's retries
kills the whole run. Consider `errorStrategy 'retry'` with a small `maxRetries` on the download
processes.

### F4. `quantify_reads_per_gene` strict validations are commented out *(low)*
The "TSV must have 2 lines / 7 columns" checks are disabled (`# Temporarily disabled`). Fine, but
"temporarily" has clearly become permanent; either delete or re-enable.

### F5. PE `--un-conc` + `cat R1 R2 > unaligned.fastq` can double-feed spike *(low)*
In PE, unaligned mates are written per-mate and concatenated into one SE file for the spike
alignment. For a pair where one mate aligned and the other didn't, behaviour depends on bowtie2's
`--un-conc` semantics; the spike count can include partial-pair reads. Minor, but it feeds siCPM.

---

## G. Numerical / edge-case details worth checking

- **Negative-strand sign convention** is threaded consistently (mirrored `-scale -1` in `06`,
  `abs()` in `08`/`10`/`11`/`16`), which is good — but it means *every* consumer must remember to
  abs(). `functional_regions.py` keeps negative values in `PROseq_3pnt.bed` col4 and only abs()es
  at summary time; any new consumer of that intermediate will get signed values.
- **`pausing_index.tsv` is 8 columns**, but `11`'s validator warns "expected 7" and the
  MEDIAN_PI awk has branch logic for both 7- and 8-column cases — leftover from a schema change.
- **Active-gene set keyed by gene *name*** (`functional_regions.py`): non-unique symbols
  (paralogs, readthroughs) collide, so one active gene can activate its namesakes' regions.
- **`min(user, auto)` for spike control**, `min()` for body offset, `max()`/`min()` sign flips —
  several of these "take the smaller/larger" choices are the kind of thing that's easy to get
  backwards; B1 is the one that actually bites.

---

## H. Documentation drift (low risk, high confusion)

- **H1.** `09_…nf` README and `detect_divergent_transcription.py` docstrings still describe the
  *old* auto-calibration (95th percentile, 10× sum, ±5 kb) while the live defaults are 65–75th
  percentile, 1.5–3× sum. The README also says divergent uses the **allMap** 3′ tracks, but
  `main.nf` actually wires the **main** 3′ tracks (`bw3p_pair_ch`). Pick one and make code+docs
  agree.
- **H2.** STEP labels in `main.nf` logging are off-by-one in places (e.g. the functional-regions
  block logs "STEP 12 | COMPLETE"; the header comments number steps 1-15 but the body goes to 17).
- **H3.** `07`'s README claims gene-level counts (see D4); `11`'s README claims 19-col output that
  matches, but the pausing schema note (7 vs 8) is stale.

---

## Suggested triage order

1. **A1** (run-on broken), **A2** (condition mis-map), **A3** (siCPM control on merge) — these
   change or empty real outputs and are small, contained fixes.
2. **C1/C2** — decide and enforce one PE mate convention, and make tracks + gene-metrics use the
   same BAM.
3. **B1** — re-examine the body-offset auto-calibration; it likely biases every pausing index.
4. **B2/B4** — make a conscious, documented decision about the divergent caller and the spike
   alignment strategy (these are *methods* choices, not bugs).
5. **E** — delete the dead modules/scripts so the next person isn't misled.
6. **F1/F2** — move report inputs and SRR caching onto Nextflow channels/`-resume` to remove the
   publish-path and string-sync fragility.

I did not modify any code (per your standing rule). Everything above is observation only; happy to
turn any single item into a concrete patch on the dev repo when you want.
