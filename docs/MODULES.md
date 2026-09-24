# Pipeline Modules

[← Back to README](../README.md)

TrackTx runs 17 modules in sequence. Here is what each one does.

**01 — Download Genome Annotations**
Downloads GTF gene annotation files from Ensembl, RefSeq, or GENCODE for the chosen reference genome. Annotations drive gene-body boundary calls, divergent transcription pairing, and functional region assignment throughout the pipeline.

**02 — Download SRA Samples**
Retrieves FASTQ files from NCBI SRA or the European Nucleotide Archive (ENA) when `sample_source: srr` is set. Uses `fasterq-dump` with multi-threaded conversion and falls back to parallel ENA downloads automatically. Raw FASTQs are stored via Nextflow's `storeDir` at `{output_dir}/00_sra_cache/{SRR}/`, so re-downloads are skipped as long as those files exist — even after the `work/` directory has been deleted between runs.

**03 — Preprocess and Quality Filter Reads**
Trims adapter sequences, removes barcodes, and applies a minimum length filter in a single `cutadapt` pass. UMI extraction with `umi_tools` follows if enabled. FastQC runs on both raw and cleaned reads to confirm quality before alignment.

**04 — Download Genome and Build Alignment Index**
Downloads the reference (and spike-in) genome FASTA, then builds a `bowtie2` alignment index. Pre-built indexes are cached and reused across runs.

**05 — Align Reads to Genome**
Aligns cleaned reads to the reference genome with `bowtie2` (end-to-end mode). For PRO-seq, R1 is reverse-complemented and paired-end data is aligned with `--ff` (`-1` original R2, `-2` RC(R1)) so the nascent-RNA 3′ end maps to the polymerase position. Spike-in reads (the genome-unaligned set) are aligned separately to derive normalization factors and are **not** run in multimapping mode. Outputs include coordinate-sorted, indexed BAM files for primary, all-mapping, and spike-in alignments.

*Multimapping (`align.multimap_k`, default 4):* bowtie2 reports up to *N* alignments per read (`-k N`). The full set becomes the **allMap** BAM — multimappers are kept at every reported locus, which makes signal visible in repetitive regions where a single best position is misleading. Filtering to one primary alignment per read (`-F 260`) gives the **main** BAM that drives counts, divergent detection, and functional-region calling. Because `-k` makes bowtie2 set MAPQ to 255, deterministic `NH:i` tags are added right after alignment (no second alignment pass), and every downstream "unique read" step selects `NH == 1`. Set `align.multimap_k: 1` for legacy single-best behaviour (then allMap ≡ main and uniqueness falls back to MAPQ).

BAMs are published to `{output_dir}/02_alignments/{sample_id}/`. To force realignment for a sample, delete its folder: `rm -rf {output_dir}/02_alignments/{sample_id}`.

**05b — Check and Merge Replicates** *(optional)*
When replicate merging is enabled (`replicates.merge: true`), performs a pairwise Pearson correlation check across BAMs using `deepTools multiBamSummary`. Replicate groups that meet the concordance threshold are merged into a single BAM before coverage track generation, with a concordance TSV written for the cohort report.

**06 — Generate Coverage Tracks**
Produces strand-specific 3′-end and 5′-end coverage tracks in bedGraph format using `bedtools genomecov`, for both the **main** BAM (best alignment per read) and the **allMap** BAM (all reported alignments). Negative-strand tracks are mirrored with `-scale -1`, and each bedGraph is sorted inline and converted to BigWig. For paired-end libraries, only Read 2 (the RC(R1) mate carrying the nascent-RNA 3′ end = polymerase position) is used for coverage, so the other mate's end cannot contaminate the tracks — this mate filtering is applied consistently to all four track sets (main/allMap × 3′/5′). For single-end data the full read set is used.

**07 — Collect Library Sizes**
Collects per-sample library sizes with `samtools idxstats`: total mapped reads in the main BAM, the allMap BAM, and the spike-in BAM. These totals form the counts master file that drives CPM and siCPM normalization-factor calculation in the next module. (Per-gene read counting is not done here — that happens in module 11 directly on the alignments.)

**08 — Normalize Coverage Tracks**
Scales raw bedGraph signal to CPM (counts per million mapped reads) and siCPM (spike-in CPM) using pre-computed factors. Positive and negative strand tracks are normalized in parallel. Outputs both bedGraph and BigWig formats for all track sets (3p, 5p, main, and allMap).

**09 — Detect Divergent Transcription**
The statistical core of the pipeline. Operates on the **main** signal track (3′ for PRO-seq, 5′ for GRO-seq). Pairs upstream antisense peaks with downstream sense peaks, computes a suite of features (signal balance, local enrichment, strand specificity), and fits a two-component Gaussian Mixture Model to separate signal from noise. Divergent regions passing the FDR threshold are written as a BED file with confidence scores. No manual thresholds are required — set `divergent_threshold: auto` and the calibration percentile handles it.

**10 — Assign Signal to Functional Regions**
Assigns normalized coverage to a hierarchical set of genomic functional regions: active promoters, gene bodies, cleavage and polyadenylation sites, enhancers, termination windows, and non-localized signal. Each position is assigned to exactly one region by sequential masking, so the categories are mutually exclusive.

**11 — Calculate Polymerase Occupancy Metrics**
Computes two complementary views of Pol II occupancy. The density metrics approach sums normalized bedGraph signal over each functional region. The gene metrics approach operates on the filtered BAM directly, computing per-gene TSS-window and gene-body coverage from which pausing indices (PI = TSS density / body density) are derived. Both approaches run in parallel so neither waits on the other. Read counting is strand-specific so only sense-strand reads contribute to each gene's TSS and body counts, eliminating contamination from antisense transcription at convergent loci, and restricted to uniquely-mapped reads (`NH == 1` when `align.multimap_k > 1`, otherwise MAPQ ≥ `pol.mapq`) so ambiguous multimappers do not inflate gene quantification. The gene-body offset is automatically calibrated from the gene-length distribution in the annotation (25th-percentile-based), so the pipeline works correctly for compact genomes such as *D. melanogaster* and *C. elegans* without manual parameter tuning.

**12 — Summarize Polymerase Metrics**
Aggregates per-sample Pol II metrics across the cohort into summary TSVs — pausing index distributions, region density tables, and normalization factor comparisons — for use in the cohort report.

**13 — Quality Control Aligned Reads**
Calculates per-sample alignment QC: total and mapped read counts, duplicate rates, unique-read rate (`NH == 1` in multimapping mode, MAPQ ≥ `qc.mapq` otherwise — reported as `uniqueness_method` in the QC JSON), strand balance (critical for PRO-seq validation), fragment length distribution (PE only), and mean genome coverage depth. Results feed the per-sample HTML reports and cohort outlier detection.

**14 — Generate Per-Sample Reports**
Produces an interactive HTML report for each sample, summarising QC metrics (including the unique-read count with its method label, and a multimapper % = 1 − unique/mapped), coverage distributions, divergent transcription statistics, and Pol II pausing results, with inline visualizations. Track links distinguish the `main` (best-alignment) and `allMap` (multimapper-aware) BigWigs so the right track is used for each purpose.

**15 — Combine Reports into Cohort**
The final step — runs after module 16 so it can incorporate signal-QC outputs into the landing page. Merges all per-sample JSON reports into a cohort-level HTML dashboard (global_summary.html) covering: by-condition QC comparisons, mapping uniqueness and multimapper % per sample, divergent transcription patterns, Pol II pausing distributions, functional region composition, normalization factor validation, replicate consistency (coefficient of variation), and an interactive sample metrics table. Also generates a modern `index.html` landing page at the output root that embeds the run-on efficiency table, KPI stats, and links to all outputs.

**16 — Cohort QC and Visualization**
Cohort-level signal QC module that runs after all per-sample tracks are ready, and before module 15 so its outputs feed the landing page. Produces: (1) a **MultiQC** HTML report aggregating all alignment logs, flagstats, and trimming statistics into a single QC dashboard; (2) **deepTools** PCA plot and Pearson correlation heatmap computed from CPM-normalized 3′ BigWigs using 10 kb genome-wide bins; (3) an **IGV session XML** file that loads all sample tracks in one click, colour-coded by condition; and (4) a **run-on efficiency table** reporting the median 5′/3′ bedGraph signal ratio across long gene bodies per sample — values close to 1.0 indicate efficient NRO run-on.
