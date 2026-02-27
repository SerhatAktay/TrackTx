// ============================================================================
// TrackTx — PRO-seq Nascent RNA Analysis Pipeline
// ============================================================================
//
// Purpose:
//   Complete PRO-seq analysis from FASTQ to functional regions and Pol-II metrics
//
// Pipeline Steps:
//   1.  Download genome annotations (GTF)
//   2.  Download SRA samples (optional)
//   3.  Preprocess and quality-filter reads
//   4.  Download genome and build alignment index (primary + spike-in)
//   5.  Align reads to genome (Bowtie2)
//   6.  Generate coverage tracks (bedGraph/BigWig)
//   7.  Quantify reads per gene
//   8.  Normalize coverage tracks (CPM/siCPM)
//   9.  Detect divergent transcription (edge-to-edge distance)
//   10. Assign signal to functional regions
//   11. Calculate polymerase occupancy metrics
//   12. Summarize polymerase metrics
//   13. Quality control aligned reads
//   14. Generate per-sample reports
//   15. Combine reports into cohort
//
// Key Updates:
//   • Divergent transcription uses edge-to-edge distance (not center-to-center)
//   • Uses primary/unique mapper tracks (not allMap) for clean detection
//   • Overlap-aware pairing algorithm (overlapping peaks always paired)
//   • Default pairing window reduced to 500bp (was 800bp)
//   • Enhanced progress monitoring and validation throughout
//
// Inputs:
//   - Samplesheet CSV (sample metadata + file paths or SRA accessions)
//   - Reference genome (hg38, mm39, or custom FASTA)
//   - Spike-in genome (optional)
//
// Outputs:
//   ${params.output_dir}/
//     ├── 01_annotations/       — GTF, genes, TSS, TES
//     ├── 02_alignments/        — BAM files, QC metrics
//     ├── 03_genome_tracks/     — Raw bedGraphs and BigWigs
//     ├── 04_counts/            — Read counts
//     ├── 05_normalized_tracks/ — CPM/siCPM tracks
//     ├── 06_divergent_tx/      — Divergent transcription calls
//     ├── 07_functional_regions/— Functional region annotations
//     ├── 08_pol_metrics/      — Per-gene polymerase occupancy metrics
//     ├── 09_pol_aggregate/    — Cohort-level summaries
//     ├── 10_qc/                — Quality control reports
//     └── 11_reports/           — HTML reports
//
// ============================================================================

nextflow.enable.dsl = 2

// ============================================================================
// PARAMETER DEFAULTS (Fallbacks for Module-Specific Parameters)
// ============================================================================
// Note: Core parameters (output_dir, samplesheet, etc.) are in nextflow.config
// These are module-specific defaults that may not be in params.yaml

// Report generation parameters
params.reports_plots = params.reports_plots ?: 0

// Conda environment overrides (if not using containers)
params.conda_norm = params.conda_norm ?: null
params.conda_pol = params.conda_pol ?: null

// QC parameters (used by preprocess_and_quality_filter_reads and quality_control_aligned_reads modules)
params.qc = params.qc instanceof Map ? params.qc : [:]
params.qc.enabled = params.qc.enabled != null ? params.qc.enabled : true
params.qc.mapq = params.qc.mapq ?: 10
params.qc.dedup = params.qc.dedup != null ? params.qc.dedup : true
params.qc.depth_max_cov = params.qc.depth_max_cov ?: 0
// Optional fail thresholds: fail process if metric below value (null = do not fail)
params.qc.fail_map_rate_below = params.qc.fail_map_rate_below  // e.g. 50 to fail if mapping rate < 50%
params.qc.fail_strand_min = params.qc.fail_strand_min           // e.g. 0.30 to fail if either strand < 30%

// Polymerase aggregate parameters (used by summarize_polymerase_metrics module)
params.pol = params.pol instanceof Map ? params.pol : [:]
params.pol.top_n = params.pol.top_n ?: 100
params.pol.plots = params.pol.plots != null ? params.pol.plots : true

// ============================================================================
// HELP MESSAGE
// ============================================================================

if (params.help) {
  log.info """
  ════════════════════════════════════════════════════════════════════════
  TrackTx Pipeline — PRO-seq Nascent RNA Analysis
  ════════════════════════════════════════════════════════════════════════
  
  Usage:
    nextflow run main.nf -entry TrackTx -profile docker \\
      --samplesheet samples.csv \\
      --reference_genome hg38 \\
      --output_dir results/
  
  Required Parameters:
    --samplesheet <CSV>           Samplesheet with sample metadata
    --reference_genome <ID>       Reference genome (hg38, mm39, or 'other')
    --output_dir <DIR>            Output directory path
  
  Optional Parameters:
    --paired_end                  Enable paired-end mode (default: false)
    --sample_source <srr|local>   Sample source (default: local)
    --spikein_genome <ID>         Spike-in genome (optional)
    --genome_fasta <FASTA>        Custom genome (if reference_genome=other)
    --spikein_fasta <FASTA>       Custom spike-in (if spikein_genome=other)
    --debug                       Enable debug logging
    --verbose                     Enable step-by-step progress logging
  
  Profiles:
    -profile docker               Use Docker containers
    -profile conda                Use Conda environments
    -profile local                Local execution
  
  Additional Flags:
    -with-report report.html      Generate execution report
    -with-timeline timeline.html  Generate timeline
    -with-dag flowchart.pdf       Generate DAG visualization
  
  New:
    • Divergent transcription uses edge-to-edge distance algorithm
    • Primary/unique mapper tracks for clean detection
    • 500bp default pairing window (optimized for promoters)
    • Overlap-aware pairing (overlapping peaks always paired)
  
  For detailed documentation, see README.md
  ════════════════════════════════════════════════════════════════════════
  """
  exit 0
}

// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

def TIMESTAMP_START = new Date().format("yyyy-MM-dd HH:mm:ss")

log.info """
════════════════════════════════════════════════════════════════════════
TrackTx Pipeline Launch
════════════════════════════════════════════════════════════════════════
Started:          ${TIMESTAMP_START}
Output directory: ${params.output_dir}
Samplesheet:      ${params.samplesheet}
Reference genome: ${params.reference_genome}
Paired-end mode:  ${params.paired_end ?: false}
Sample source:    ${params.sample_source ?: 'local'}
Spike-in genome:  ${params.spikein_genome ?: 'None'}
Debug mode:       ${params.debug ?: false}
════════════════════════════════════════════════════════════════════════
""".stripIndent()

// Critical parameter validation
if (!params.output_dir) {
  error "PIPELINE | ERROR | Missing required parameter: --output_dir"
}
if (!params.samplesheet) {
  error "PIPELINE | ERROR | Missing required parameter: --samplesheet"
}
if (!params.reference_genome) {
  error "PIPELINE | ERROR | Missing required parameter: --reference_genome"
}

// Custom genome validation
if (params.reference_genome == 'other' && !params.genome_fasta) {
  error "PIPELINE | ERROR | When reference_genome=other, must provide --genome_fasta"
}
if (params.spikein_genome == 'other' && !params.spikein_fasta) {
  error "PIPELINE | ERROR | When spikein_genome=other, must provide --spikein_fasta"
}

// Validate samplesheet exists (resolve relative paths from projectDir)
def samplesheetPath = params.samplesheet?.trim()
def samplesheetFile = samplesheetPath && new File(samplesheetPath).isAbsolute() 
  ? file(samplesheetPath) 
  : file("${projectDir}/${samplesheetPath ?: 'samplesheet.csv'}")
if (!samplesheetFile.exists()) {
  error "PIPELINE | ERROR | Samplesheet not found: ${params.samplesheet}"
}

// Validate samplesheet has at least one data row (fail before workflow starts)
def samplesheetLines = samplesheetFile.readLines()
def dataRows = samplesheetLines.drop(1).findAll { it.trim() }
if (dataRows.isEmpty()) {
  error "PIPELINE | ERROR | Samplesheet has no data rows. Expected format: sample,condition,timepoint,replicate,file1,file2"
}

// Fail-fast: validate input files exist for local samples (prevents silent skip → confusing Step 13 error)
if (params.sample_source != 'srr') {
  def header = samplesheetLines[0]?.split(',')?.collect { it.trim() }
  def sampleIdx = header?.findIndexOf { it?.toLowerCase() == 'sample' }
  def file1Idx = header?.findIndexOf { it?.toLowerCase() == 'file1' }
  def file2Idx = header?.findIndexOf { it?.toLowerCase() == 'file2' }
  if (sampleIdx != null && sampleIdx >= 0 && file1Idx != null && file1Idx >= 0) {
    def missing = []
    dataRows.eachWithIndex { line, i ->
      def cols = line.split(',', -1).collect { it?.trim() }
      def sample = cols.size() > sampleIdx ? cols[sampleIdx] : ''
      def resolve = { String p ->
        if (!p?.trim()) return null
        def s = p.trim()
        def f = new File(s)
        return f.isAbsolute() ? f : new File(projectDir.toString(), s)
      }
      def f1 = cols.size() > file1Idx ? resolve(cols[file1Idx]) : null
      def f2 = (params.paired_end && file2Idx != null && file2Idx >= 0 && cols.size() > file2Idx) ? resolve(cols[file2Idx]) : null
      if (f1 && !f1.exists()) missing << "${sample}: file1 not found: ${f1}"
      if (f2 && !f2.exists()) missing << "${sample}: file2 not found: ${f2}"
    }
    if (!missing.isEmpty()) {
      error """PIPELINE | ERROR | Input files not found (${missing.size()} issue(s)).
  Fix paths in samplesheet or ensure files exist before running.

  Missing/not found:
  ${missing.take(10).join('\n  ')}${missing.size() > 10 ? '\n  ... and ' + (missing.size() - 10) + ' more' : ''}

  Expected format: sample,condition,timepoint,replicate,file1,file2
  Paths are relative to: ${projectDir}"""
    }
  }
}

if (params.verbose) log.info "PIPELINE | VALIDATE | Parameter validation complete"

// ============================================================================
// MODULE IMPORTS
// ============================================================================

def MOD = "${projectDir}/modules"

include { download_genome_annotations                             } from "${MOD}/01_download_genome_annotations.nf"
include { download_sra_samples                                   } from "${MOD}/02_download_sra_samples.nf"
include { preprocess_and_quality_filter_reads                    } from "${MOD}/03_preprocess_and_quality_filter_reads.nf"
include { download_genome_and_build_alignment_index as build_index } from "${MOD}/04_download_genome_and_build_alignment_index.nf"
include { download_genome_and_build_alignment_index as spike_index } from "${MOD}/04_download_genome_and_build_alignment_index.nf"
include { align_reads_to_genome                                  } from "${MOD}/05_align_reads_to_genome.nf"
include { generate_coverage_tracks                               } from "${MOD}/06_generate_coverage_tracks.nf"
include { quantify_reads_per_gene                                } from "${MOD}/07_quantify_reads_per_gene.nf"
include { normalize_coverage_tracks                              } from "${MOD}/08_normalize_coverage_tracks.nf"
include { detect_divergent_transcription                         } from "${MOD}/09_detect_divergent_transcription.nf"
include { assign_signal_to_functional_regions                    } from "${MOD}/10_assign_signal_to_functional_regions.nf"
include { calculate_polymerase_occupancy_metrics                 } from "${MOD}/11_calculate_polymerase_occupancy_metrics.nf"
include { summarize_polymerase_metrics                           } from "${MOD}/12_summarize_polymerase_metrics.nf"
include { quality_control_aligned_reads                          } from "${MOD}/13_quality_control_aligned_reads.nf"
include { generate_per_sample_reports                            } from "${MOD}/14_generate_per_sample_reports.nf"
include { combine_reports_into_cohort                            } from "${MOD}/15_combine_reports_into_cohort.nf"

if (params.verbose) log.info "PIPELINE | IMPORT | All modules loaded successfully"

// Resolve local FASTQ path: use as-is if absolute, else relative to projectDir
def resolveLocalPath(path) {
  if (!path?.trim()) return null
  def p = path.trim()
  def f = new File(p)
  return f.isAbsolute() ? file(p) : file("${projectDir}/${p}")
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow TrackTx {

  if (params.verbose) {
    log.info "════════════════════════════════════════════════════════════════════════"
    log.info "WORKFLOW | Starting TrackTx Analysis"
    log.info "════════════════════════════════════════════════════════════════════════"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 1: Download Annotations
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 1 | Download Genome Annotations"
    log.info "─".multiply(80)
  }
  
  download_genome_annotations()
  
  gtf_ch   = download_genome_annotations.out[0]
  genes_ch = download_genome_annotations.out.genes
  tss_ch   = download_genome_annotations.out.tss
  tes_ch   = download_genome_annotations.out.tes
  


  // ══════════════════════════════════════════════════════════════════════════
  // STEP 2: Parse Samplesheet
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 2 | Parse Samplesheet"
    log.info "─".multiply(80)
    log.info "STEP 2 | INPUT | File: ${params.samplesheet}"
  }

  // Use .tap() to duplicate channel for counting without consuming it
  samples_ch = Channel
    .fromPath(samplesheetFile)
    .splitCsv(header: true)
    .map { row ->
      // Validate required fields
      if (!row.sample || !row.file1) {
        error "STEP 2 | ERROR | Samplesheet row missing 'sample' or 'file1': ${row}"
      }
      
      // Validate file2 for paired-end local FASTQ mode
      if (params.sample_source != 'srr' && params.paired_end && !row.file2) {
        error "STEP 2 | ERROR | Paired-end mode requires 'file2' for sample ${row.sample}"
      }

      // Build sample ID with replicate
      def sample_id = "${row.sample}_r${row.replicate ?: 1}"
      
      // Build reads list based on source (supports absolute or relative paths for local files)
      def reads = (params.sample_source == 'srr')
        ? [row.file1.trim(), row.file2?.trim()]
        : (params.paired_end 
            ? [resolveLocalPath(row.file1), resolveLocalPath(row.file2)] 
            : [resolveLocalPath(row.file1)])
      
      tuple(
        sample_id,
        reads,
        (row.treatment ?: row.condition ?: row.sample) ?: '',  // condition (treatment/condition column)
        row.timepoint ?: '',
        row.replicate ?: 1
      )
    }
    .tap { samples_for_count }

  // Count samples without consuming main channel; fail fast if none parsed
  samples_for_count.count().subscribe { count ->
    if (count == 0) {
      error "PIPELINE | ERROR | No samples parsed from samplesheet. Check column names (sample, file1), delimiter (comma), and encoding (UTF-8)."
    }
    if (params.verbose) log.info "STEP 2 | COMPLETE | Loaded ${count} samples from samplesheet"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 3: Download SRR Data (optional)
  // ══════════════════════════════════════════════════════════════════════════
  
  prepared_input_ch = null
  
  if (params.sample_source == 'srr') {
    if (params.verbose) {
      log.info "─".multiply(80)
      log.info "STEP 3 | Download SRR Data from NCBI"
      log.info "─".multiply(80)
    }
    
    // Explicitly evaluate parameter to avoid closure comparison issues
    def isPairedEnd = params.paired_end ? true : false
    
    download_sra_samples(
      samples_ch.map { sid, reads, c, t, r ->
        tuple(sid, reads[0], c, t, r, isPairedEnd)
      }
    )
    
    prepared_input_ch = download_sra_samples.out[0].map { sid, fq1, fq2, c, t, r ->
      tuple(sid, [file(fq1), file(fq2)].findAll(), c, t, r)
    }
    
    // Monitor downloads
    prepared_input_ch.subscribe { sid, _reads, _cond, _time, _rep ->
      if (params.verbose) log.info "STEP 3 | DOWNLOAD | ${sid} complete"
    }
    

    
  } else {
    if (params.verbose) {
      log.info "─".multiply(80)
      log.info "STEP 3 | Using Local FASTQ Files"
      log.info "─".multiply(80)
    }
    prepared_input_ch = samples_ch

  }
  
  // ══════════════════════════════════════════════════════════════════════════
  // STEP 4: Prepare and Clean Input Reads
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 4 | Preprocess Input Reads"
    log.info "─".multiply(80)
    log.info "STEP 4 | CONFIG | Mode: ${params.paired_end ? 'Paired-end' : 'Single-end'}"
    log.info "STEP 4 | CONFIG | Trimming and quality filtering enabled"
  }
  
  (clean_fastq_ch, fastqc_ch) = preprocess_and_quality_filter_reads(
    prepared_input_ch,
    Channel.value(params.paired_end ? 'PE' : 'SE')
  )

  // Monitor preprocessing - FIX: _r used twice before!
  clean_fastq_ch.subscribe { sid, _r1, _r2opt, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 4 | CLEAN | ${sid} preprocessing complete"
  }

  // Create sentinel files for SE mode (R2 placeholder)
  def assetsDir = params.assets_dir ?: "${projectDir}/assets"
  new File(assetsDir).mkdirs()
  
  def noR2Path = "${assetsDir}/NO_R2.fastq"
  def noBGPath = "${assetsDir}/EMPTY.bedgraph"
  def noBGPosPath = "${assetsDir}/EMPTY_POS.bedgraph"
  def noBGNegPath = "${assetsDir}/EMPTY_NEG.bedgraph"
  // Unique placeholders for normalize_coverage_tracks to avoid staging collisions
  def noBG5pPosPath = "${assetsDir}/EMPTY_5P_POS.bedgraph"
  def noBG5pNegPath = "${assetsDir}/EMPTY_5P_NEG.bedgraph"
  def noBGAm5pPosPath = "${assetsDir}/EMPTY_AM5P_POS.bedgraph"
  def noBGAm5pNegPath = "${assetsDir}/EMPTY_AM5P_NEG.bedgraph"
  // Placeholder for align when no spike-in (avoids Channel.empty() causing 0 invocations)
  def noSpikeFaPath = "${assetsDir}/EMPTY_SPIKE.fa"
  
  // Create all placeholder files if they don't exist
  [noR2Path, noBGPath, noBGPosPath, noBGNegPath, 
   noBG5pPosPath, noBG5pNegPath, noBGAm5pPosPath, noBGAm5pNegPath].each { path ->
    if (!new File(path).exists()) {
      new File(path).text = ''
    }
  }
  if (!new File(noSpikeFaPath).exists()) {
    new File(noSpikeFaPath).text = ">none\nN\n"
  }

  // Ensure R2 is always present (use sentinel for SE)
  clean_fastq_with_r2 = clean_fastq_ch.map { sid, r1, r2opt, c, t, r ->
    def r2_file = (r2opt && file(r2opt).exists() && file(r2opt).size() > 0) 
      ? file(r2opt) 
      : file(noR2Path)
    tuple(sid, r1, r2_file, c, t, r)
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 5: Build Genome Indices
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 5 | Build Bowtie2 Indices"
    log.info "─".multiply(80)
  }

  // Primary genome
  reference_fa = params.reference_genome == 'other'
    ? file(params.genome_fasta)
    : file("${projectDir}/genomes/${params.reference_genome}.fa")

  if (params.verbose) log.info "STEP 5 | INDEX | Primary genome: ${params.reference_genome}"

  build_index(
    Channel.value(tuple(
      params.reference_genome == 'other' ? 'custom' : params.reference_genome,
      params.reference_genome == 'other' ? 'custom' : 'ucsc',
      reference_fa
    ))
  )
  
  ref_meta_ch = build_index.out.ref_meta
  ref_idx_ch  = build_index.out.index_files

  if (params.verbose) log.info "STEP 5 | INDEX | Primary genome index built"

  // Spike-in genome (optional)
  // CRITICAL: When no spike-in, use placeholder channels (not Channel.empty()) so align_reads
  // receives one value per input and runs. Empty channels cause 0 invocations → downstream fails at Step 13.
  def noSpikeFa = file(noSpikeFaPath)
  spike_meta_ch = Channel.value(tuple('none', 'none', noSpikeFa))
  spike_idx_ch  = Channel.value(noSpikeFa)
  
  if (params.spikein_genome && params.spikein_genome != 'None') {
    if (params.verbose) log.info "STEP 5 | INDEX | Spike-in genome: ${params.spikein_genome}"
    
    spike_fa = params.spikein_genome == 'other'
      ? file(params.spikein_fasta)
      : file("${projectDir}/genomes/${params.spikein_genome}.fa")

    spike_index(
      Channel.value(tuple(
        params.spikein_genome == 'other' ? 'custom_spike' : params.spikein_genome,
        params.spikein_genome == 'other' ? 'custom' : 'ucsc',
        spike_fa
      ))
    )
    
    spike_meta_ch = spike_index.out.ref_meta
    spike_idx_ch  = spike_index.out.index_files
    
    if (params.verbose) log.info "STEP 5 | INDEX | Spike-in index built"
  } else {
    if (params.verbose) log.info "STEP 5 | INDEX | No spike-in genome specified (using placeholder)"
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 6: Alignment
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 6 | Align Reads with Bowtie2"
    log.info "─".multiply(80)
    log.info "STEP 6 | CONFIG | Aligner: Bowtie2"
    log.info "STEP 6 | CONFIG | Mode: ${params.paired_end ? 'Paired-end' : 'Single-end'}"
  }

  align_reads_to_genome(
    clean_fastq_with_r2,
    ref_meta_ch,
    ref_idx_ch,
    spike_meta_ch,
    spike_idx_ch,
    params.paired_end ?: false
  )
  
  aligned_ch = align_reads_to_genome.out[0]
  
  // Monitor alignments
  aligned_ch.subscribe { sid, _bam, _allbam, _spike, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 6 | ALIGN | ${sid} alignment complete"
  }
  
  // aligned_ch: (sample_id, sample.bam, sample_allMap.bam, spikein.bam, condition, timepoint, replicate)



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 7: Generate Coverage Tracks
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 7 | Generate Coverage Tracks"
    log.info "─".multiply(80)
    log.info "STEP 7 | CONFIG | Formats: bedGraph, BigWig"
    log.info "STEP 7 | CONFIG | Orientations: 3' end (primary), 5' end (optional)"
    log.info "STEP 7 | CONFIG | Strands: Positive, Negative (separate)"
  }

  genome_fa_ch = ref_meta_ch.map { id, prov, fa -> fa }

  tracks_input_ch = aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
    tuple(sid, filt_bam, spike_bam, c, t, r)
  }
  
  allmap_bam_ch = aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
    all_bam
  }

  generate_coverage_tracks(
    tracks_input_ch,
    genome_fa_ch,
    Channel.value(params.paired_end ? 'true' : 'false'),
    allmap_bam_ch
  )

  bw3p_pair_ch    = generate_coverage_tracks.out.bw3p_pair      // PRIMARY/UNIQUE MAPPERS ✓
  bw5p_pair_ch    = generate_coverage_tracks.out.bw5p_pair
  allmap3p_pair_ch = generate_coverage_tracks.out.allmap3p_pair  // MULTIMAPPERS (not for divergent tx)
  allmap5p_pair_ch = generate_coverage_tracks.out.allmap5p_pair
  tracks_ch       = generate_coverage_tracks.out.track_tuple
  dedup_stats_ch  = generate_coverage_tracks.out.dedup_stats

  // Monitor track generation
  tracks_ch.subscribe { sid, _fb, _sb, p3, n3, _p5, _n5, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 7 | TRACKS | ${sid} → 3' tracks generated"
  }

  if (params.debug) {
    tracks_ch.view { sid, fb, sb, p3, n3, p5, n5, c, t, r ->
      "DEBUG | TRACKS | ${sid} | 3p: ${p3}, ${n3} | 5p: ${p5 ?: 'N/A'}, ${n5 ?: 'N/A'}"
    }
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 8: Collect Read Counts
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 8 | Collect Read Counts"
    log.info "─".multiply(80)
    log.info "STEP 8 | PURPOSE | Read counts for CPM/siCPM normalization"
  }

  counts_tsvs = quantify_reads_per_gene(
    aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
      tuple(sid, filt_bam, all_bam, spike_bam ?: '-', c, t, r)
    }
  )

  counts_master = counts_tsvs.counts
    .map { sid, tsv, c, t, r -> tsv }
    .collectFile(
      name: 'counts_master.tsv',
      storeDir: "${params.output_dir}/04_counts",
      keepHeader: true,
      newLine: true
    )

  counts_master.subscribe {
    if (params.verbose) log.info "STEP 8 | COMPLETE | Master counts file: ${it}"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 9: Normalize Tracks
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 9 | Normalize Tracks"
    log.info "─".multiply(80)
    log.info "STEP 9 | CONFIG | Methods: CPM (counts per million), siCPM (spike-in CPM)"
    log.info "STEP 9 | CONFIG | Tracks: 3' primary, 3' allMap"
  }

  // Prepare inputs - remove redundant metadata from joins for clarity
  // NOTE: Must wrap paths with file() to ensure proper staging when using path() inputs
  // NOTE: Use unique placeholder files to avoid Nextflow staging name collisions
  norm_main_kv = tracks_ch.map { sid, bam, spk, p3, n3, p5, n5, c, t, r ->
    def p3_file = file(p3)
    def n3_file = file(n3)
    def p5_file = (p5 && file(p5).exists() && file(p5).size() > 0) ? file(p5) : file(noBG5pPosPath)
    def n5_file = (n5 && file(n5).exists() && file(n5).size() > 0) ? file(n5) : file(noBG5pNegPath)
    tuple(sid, tuple(p3_file, n3_file, p5_file, n5_file, c, t, r))
  }

  norm_allmap3_kv = allmap3p_pair_ch.map { sid, ap3, an3, bwp, bwn, c, t, r ->
    tuple(sid, tuple(file(ap3), file(an3)))
  }

  norm_allmap5_kv = allmap5p_pair_ch.map { sid, ap5, an5, bwp, bwn, c, t, r ->
    def ap5_file = (ap5 && file(ap5).exists() && file(ap5).size() > 0) ? file(ap5) : file(noBGAm5pPosPath)
    def an5_file = (an5 && file(an5).exists() && file(an5).size() > 0) ? file(an5) : file(noBGAm5pNegPath)
    tuple(sid, tuple(ap5_file, an5_file))
  }

  // Join using sample_id as key only, rebuild full tuple at end
  norm_input_ch = norm_main_kv
    .join(norm_allmap3_kv)
    .join(norm_allmap5_kv)
    .combine(counts_master)
    .combine(genes_ch)
    .map { sid, main_data, allmap3_data, allmap5_data, cm, genes ->
      def (p3, n3, p5, n5, c, t, r) = main_data
      def (ap3, an3) = allmap3_data
      def (ap5, an5) = allmap5_data
      tuple(sid, p3, n3, p5, n5, ap3, an3, ap5, an5, c, t, r, cm, genes)
    }

  normalize_coverage_tracks(norm_input_ch, genome_fa_ch)
  
  norm_tracks_ch  = normalize_coverage_tracks.out.norm_tuple
  norm_factors_ch = norm_tracks_ch.map { sid, p3, n3, nf, c, t, r ->
    tuple(sid, nf, c, t, r)
  }

  // Monitor normalization
  norm_tracks_ch.subscribe { sid, _p3, _n3, _nf, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 9 | NORMALIZE | ${sid} normalization complete"
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 10: Detect Divergent Transcription
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "═".multiply(80)
    log.info "STEP 10 | Detect Divergent Transcription"
    log.info "═".multiply(80)
    log.info "STEP 10 | CONFIG | Algorithm: Edge-to-edge distance with overlap awareness"
    log.info "STEP 10 | CONFIG | Default window: 500bp (edge-to-edge gap, not center-to-center)"
    log.info "STEP 10 | CONFIG | Track source: PRIMARY/UNIQUE mappers (filtered BAM)"
    log.info "STEP 10 | CONFIG | NOT using allMap (multimapper) tracks"
    log.info "STEP 10 | CONFIG | Pairing logic: Overlapping peaks always paired, then gap ≤ window"
  }

  // CRITICAL FIX: Use filtered/primary tracks (unique mappers), NOT allMap
  // 
  // bw3p_pair_ch     → from sample.bam (filtered = primary/unique mappers) ✓ CORRECT
  // allmap3p_pair_ch → from sample_allMap.bam (includes multimappers)     ✗ WRONG
  //
  // For clean divergent transcription detection, we want primary/unique mappers
  // to avoid inflated signal in repetitive regions from multimapping reads.
  //
  // CACHE FIX: Simple map without groupTuple or file checks to avoid cache invalidation
  divergent_input_ch = bw3p_pair_ch  // ✓ PRIMARY/UNIQUE MAPPERS
    .map { sid, pos3_bg, neg3_bg, bwp, bwn, c, t, r ->
      tuple(sid, pos3_bg, neg3_bg, c, t, r)
    }
  
  // Add progress monitoring
  divergent_input_ch.subscribe { sid, _pos, _neg, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 10 | INPUT | Sample ready for detection: ${sid}"
  }

  // Pass parameters explicitly for better cache control
  detect_divergent_transcription(
    divergent_input_ch,
    Channel.value(params.advanced?.divergent_threshold ?: 'auto'),
    Channel.value(params.advanced?.divergent_sum_thr ?: 'auto'),
    Channel.value(params.advanced?.divergent_fdr ?: 0.08),
    Channel.value(params.advanced?.divergent_nt_window ?: 1000),
    Channel.value(params.advanced?.divergent_balance ?: 0.0),
    Channel.value(params.advanced?.divergent_bin_gap ?: 100),
    Channel.value(params.advanced?.divergent_calibration_percentile ?: 65.0),
    Channel.value(params.advanced?.divergent_calibration_sum_multiplier ?: 1.5),
    Channel.value(params.advanced?.divergent_calibration_background_lower ?: false),
    Channel.value(params.advanced?.divergent_merge_gap ?: 150)
  )
  
  divergent_tx_ch = detect_divergent_transcription.out.bed
  
  // Monitor completion with region counts - FIX: Three duplicate _ before!
  divergent_tx_ch.subscribe { sid, bed, _cond, _time, _rep ->
    def count = 0
    if (file(bed).exists() && file(bed).size() > 0) {
      count = file(bed).readLines().findAll { !it.startsWith('#') }.size()
    }
    if (params.verbose) log.info "STEP 10 | COMPLETE | ${sid} → ${count} divergent regions detected"
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 11: Call Functional Regions
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 11 | Call Functional Regions"
    log.info "─".multiply(80)
    log.info "STEP 11 | CONFIG | Input: Divergent transcription regions + RAW tracks"
    log.info "STEP 11 | CONFIG | Output: Promoter, enhancer, gene body annotations"
  }

  // IMPORTANT: Use RAW bedGraphs (not normalized) as per original pipeline design
  func_input_ch = divergent_tx_ch
    .map { sid, div_bed, c, t, r -> 
      tuple(sid, tuple(div_bed, c, t, r))
    }
    .join(
      tracks_ch.map { sid, bam, spk, pos3_raw, neg3_raw, p5, n5, c, t, r ->
        tuple(sid, tuple(pos3_raw, neg3_raw, c, t, r))
      }
    )
    .map { sid, div_data, track_data ->
      def (div_bed, c, t, r) = div_data
      def (pos3_raw, neg3_raw, c2, t2, r2) = track_data
      // Empty placeholders for siCPM (not used in functional regions calling)
      // FIX: Use distinct files to avoid input name collision in process
      tuple(sid, div_bed, pos3_raw, neg3_raw, file(noBGPosPath), file(noBGNegPath), c, t, r)
    }

  assign_signal_to_functional_regions(
    func_input_ch,
    gtf_ch,
    Channel.value(file("${projectDir}/bin/functional_regions.py")),
    genes_ch,
    tss_ch,
    tes_ch
  )

  functional_regions_ch = assign_signal_to_functional_regions.out.main
  
  functional_regions_bed_ch = functional_regions_ch.map { sid, bed, fsum, c, t, r ->
    tuple(sid, bed)
  }
  
  functional_regions_sum_ch = functional_regions_ch.map { sid, bed, fsum, c, t, r ->
    tuple(sid, fsum, c, t, r)
  }

  // Monitor functional region calling
  functional_regions_ch.subscribe { sid, bed, _fsum, _cond, _time, _rep ->
    def count = 0
    if (file(bed).exists() && file(bed).size() > 0) {
      count = file(bed).readLines().findAll { !it.startsWith('#') }.size()
    }
    if (params.verbose) log.info "STEP 11 | COMPLETE | ${sid} → ${count} functional regions annotated"
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 12: Calculate Pol-II Metrics
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 12 | Calculate Pol-II Metrics"
    log.info "─".multiply(80)
    log.info "STEP 12 | CONFIG | Metrics: Density, pausing index, traveling ratio"
    log.info "STEP 12 | CONFIG | Input: Normalized CPM/siCPM tracks + functional regions"
  }

  // Prepare Pol-II inputs (uses normalized CPM/siCPM tracks)
  // Use BAM from generate_coverage_tracks (same BAM as density/tracks; deduped when UMI on)
  pol_input_ch = generate_coverage_tracks.out.bam_for_tracks
    .map { sid, bam, c, t, r ->
      tuple(sid, tuple(bam, c, t, r))
    }
    .join(functional_regions_bed_ch)
    .join(
      norm_tracks_ch.map { sid, pos3_cpm, neg3_cpm, factors, c, t, r ->
        // Construct siCPM paths based on directory structure from normalize_coverage_tracks
        def normDir = "${params.output_dir}/05_normalized_tracks/${sid}"
        def pos3_sicpm = file("${normDir}/3p/${sid}.3p.pos.sicpm.bedgraph")
        def neg3_sicpm = file("${normDir}/3p/${sid}.3p.neg.sicpm.bedgraph")
        tuple(sid, tuple(pos3_cpm, neg3_cpm, pos3_sicpm, neg3_sicpm))
      }
    )
    .map { sid, bam_data, bed, norm_data ->
      def (bam, c, t, r) = bam_data
      def (pos_cpm, neg_cpm, pos_si, neg_si) = norm_data
      tuple(sid, bam, bed, pos_cpm, neg_cpm, pos_si, neg_si, c, t, r)
    }

  calculate_polymerase_occupancy_metrics(pol_input_ch, gtf_ch)
  
  pol_gene_ch    = calculate_polymerase_occupancy_metrics.out.genes
  pol_density_ch = calculate_polymerase_occupancy_metrics.out.density
  pol_pausing_ch = calculate_polymerase_occupancy_metrics.out.pausing

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 13: Summarize Pol-II Metrics (Cohort-Level)
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 13 | Summarize Pol-II Metrics"
    log.info "─".multiply(80)
    log.info "STEP 13 | PURPOSE | Cohort-level aggregation and visualization"
  }

  // Sort once and reuse to avoid channel exhaustion
  pol_sorted = pol_gene_ch
    .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by sample_id once

  // Build samples TSV - use sorted list for deterministic indexing (cache-friendly)
  samples_lines = pol_sorted
    .flatMap { sorted_list ->
      sorted_list.withIndex().collect { item, idx ->
        def (sid, genes, c, t, r) = item
        // No need to check file existence - Nextflow guarantees output files exist
        "${sid}\t${c ?: 'NA'}\t${t ?: 'NA'}\t${r ?: '1'}\tmetric_${idx + 1}"
      }
    }

  samples_tsv = Channel
    .of('sample_id\tcondition\ttimepoint\treplicate\tfile')
    .concat(samples_lines)
    .collectFile(
      name: 'samples.tsv',
      newLine: true
    )
    .map { tsv_file ->
      // Sanitize TSV (remove BOM, carriage returns, duplicate headers)
      def clean_text = tsv_file.text
        .replace('\uFEFF', '')
        .replace('\r', '')
      
      def lines = clean_text.readLines()
      def header = 'sample_id\tcondition\ttimepoint\treplicate\tfile'
      def clean_lines = [header]
      def rejected = []
      
      lines.drop(1).each { line ->
        if (line.trim() && line != header) {
          def cols = line.split('\t', -1)
          if (cols.size() == 5 && cols[4] && cols[4] != 'file') {
            clean_lines << line
          } else {
            rejected << "cols=${cols.size()} file_col='${cols.size() > 4 ? cols[4] : 'N/A'}' | ${line.take(80)}..."
          }
        }
      }
      
      if (clean_lines.size() == 1) {
        def rawCount = lines.size() - 1
        def rejectMsg = rejected ? "\n  Rejected rows (first 5): ${rejected.take(5).join('\n  ')}" : ''
        error """STEP 13 | ERROR | No valid samples in Pol-II aggregate TSV
  Raw data rows: ${rawCount} | Valid after parse: ${clean_lines.size() - 1}${rejectMsg}
  This usually means no samples reached calculate_polymerase_occupancy_metrics. Possible causes:
  - An upstream process failed (check preprocess_and_quality_filter_reads, align_reads_to_genome, generate_coverage_tracks, normalize_coverage_tracks)
  - Sample ID mismatch in join operations (pol_input joins bam_for_tracks + functional_regions + norm_tracks)
  - Pipeline was interrupted before any sample completed calculate_polymerase_occupancy_metrics
  - Run with -resume to continue from last successful run, or -resume false to rule out stale cache
  Check .nextflow.log and work/*/ for failed tasks."""
      }
      
      tsv_file.text = clean_lines.join('\n') + '\n'
      if (params.verbose) log.info "STEP 13 | INPUT | Samples TSV prepared: ${clean_lines.size() - 1} samples"
      tsv_file
    }

  // Collect pol files in same sorted order as TSV for matching indices
  pol_files_ch = pol_gene_ch
    .toSortedList { a, b -> a[0] <=> b[0] }
    .flatMap { sorted_list ->
      sorted_list.collect { sid, genes, c, t, r -> file(genes) }
    }
  
  // Combine TSV with collected files and pass to process
  summarize_polymerase_metrics(samples_tsv, pol_files_ch.collect())



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 14: Quality Control
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 14 | Quality Control Analysis"
    log.info "─".multiply(80)
    log.info "STEP 14 | CONFIG | Metrics: Mapping rates, strand bias, coverage depth"
  }

  qc_input_ch = aligned_ch
    .map { sid, bam, all, spike, c, t, r ->
      tuple(sid, tuple(bam, c, t, r))
    }
    .join(dedup_stats_ch.map { sid, dedup_stats, c, t, r ->
      tuple(sid, dedup_stats)
    })
    .map { sid, bam_data, dedup_stats ->
      def (bam, c, t, r) = bam_data
      tuple(sid, bam, dedup_stats, c, t, r)
    }

  quality_control_aligned_reads(qc_input_ch)
  
  qc_json_meta_ch = quality_control_aligned_reads.out.json_meta

  // Monitor QC completion
  qc_json_meta_ch.subscribe { sid, _json, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 14 | COMPLETE | ${sid} QC analysis finished"
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 15: Generate Per-Sample Reports
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 15 | Generate Per-Sample Reports"
    log.info "─".multiply(80)
    log.info "STEP 15 | CONFIG | Format: HTML with embedded plots"
  }

  // Helper function to resolve file paths safely
  def resolvePath = { String path ->
    def f = file(path)
    f.exists() ? f.toString() : ''
  }

  // Build comprehensive report input channel
  report_input_ch = aligned_ch
    .map { sid, bam, all, spike, c, t, r ->
      tuple(sid, tuple(c ?: 'NA', t ?: 'NA', r ?: '1'))
    }
    .join(divergent_tx_ch.map { sid, bed, c, t, r -> tuple(sid, bed) })
    .join(functional_regions_sum_ch.map { sid, fsum, c, t, r -> tuple(sid, fsum) })
    .join(pol_density_ch.map { sid, dens, c, t, r -> tuple(sid, dens) })
    .join(pol_pausing_ch.map { sid, paus, c, t, r -> tuple(sid, paus) })
    .join(norm_factors_ch.map { sid, nf, c, t, r -> tuple(sid, nf) })
    .join(dedup_stats_ch.map { sid, dedup, c, t, r -> tuple(sid, dedup) })
    .join(qc_json_meta_ch.map { sid, qc, c, t, r -> tuple(sid, qc) })
    .map { sid, meta, div_bed, fsum, dens, paus, norm, dedup, qc ->
      def (c, t, r) = meta
      
      // Resolve BigWig and bedGraph paths
      def normDir = "${params.output_dir}/05_normalized_tracks/${sid}"
      def tracksDir = "${params.output_dir}/03_genome_tracks/${sid}"
      
      def bw_pos3 = resolvePath("${normDir}/3p/${sid}.3p.pos.cpm.bw")
      def bw_neg3 = resolvePath("${normDir}/3p/${sid}.3p.neg.cpm.bw")
      def bw_allmap_pos3 = resolvePath("${normDir}/3p/${sid}.allMap.3p.pos.cpm.bw")
      def bw_allmap_neg3 = resolvePath("${normDir}/3p/${sid}.allMap.3p.neg.cpm.bw")
      
      def raw_allmap_pos3 = resolvePath("${tracksDir}/3p/${sid}.allMap.3p.pos.bedgraph")
      def raw_allmap_neg3 = resolvePath("${tracksDir}/3p/${sid}.allMap.3p.neg.bedgraph")
      
      tuple(
        sid, 
        div_bed, fsum, dens, paus, 
        norm, dedup, qc,
        raw_allmap_pos3, raw_allmap_neg3,
        bw_pos3, bw_neg3, bw_allmap_pos3, bw_allmap_neg3,
        c, t, r
      )
    }

  if (params.debug) {
    report_input_ch.view { sid, div, fsum, dens, paus, norm, dedup, qc, 
                                ap, an, bwp, bwn, bwap, bwan, c, t, r ->
      "DEBUG | REPORT | ${sid} | Files: div=${div.name}, fsum=${fsum.name}, " +
      "dens=${dens.name}, paus=${paus.name}"
    }
  }

  generate_per_sample_reports(report_input_ch)
  
  report_json_ch = generate_per_sample_reports.out[2]

  // Monitor report generation
  report_json_ch.subscribe { json_file ->
    if (params.verbose) log.info "STEP 15 | COMPLETE | Report generated: ${json_file.name}"
  }



  // ══════════════════════════════════════════════════════════════════════════
  // STEP 16: Combine Reports
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "─".multiply(80)
    log.info "STEP 16 | Combine Reports"
    log.info "─".multiply(80)
    log.info "STEP 16 | PURPOSE | Generate cohort-level summary report"
  }

  // Collect reports and stage with sequential names to avoid collisions
  per_sample_reports = report_json_ch
    .toSortedList { a, b -> a.name <=> b.name }
  
  combine_reports_into_cohort(per_sample_reports)



  // ══════════════════════════════════════════════════════════════════════════
  // PIPELINE COMPLETE
  // ══════════════════════════════════════════════════════════════════════════
  
  def TIMESTAMP_END = new Date().format("yyyy-MM-dd HH:mm:ss")
  
  if (params.verbose) {
    log.info "═".multiply(80)
    log.info "PIPELINE COMPLETE | TrackTx Analysis Finished Successfully"
    log.info "═".multiply(80)
    log.info "Started:   ${TIMESTAMP_START}"
    log.info "Finished:  ${TIMESTAMP_END}"
    log.info "Results:   ${params.output_dir}"
    log.info "═".multiply(80)
    log.info ""
    log.info "Next Steps:"
    log.info "  • Review QC reports in: ${params.output_dir}/10_qc/"
    log.info "  • Check HTML reports in: ${params.output_dir}/11_reports/"
    log.info "  • Examine divergent transcription in: ${params.output_dir}/06_divergent_tx/"
    log.info "  • Analyze Pol-II metrics in: ${params.output_dir}/08_pol_metrics/"
    log.info "═".multiply(80)
  }
}