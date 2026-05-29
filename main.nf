// ============================================================================
// TrackTx — PRO-seq Nascent RNA Analysis Pipeline
// ============================================================================
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
// ============================================================================

// ============================================================================
// PARAMETER DEFAULTS (top-level param declarations — allowed by strict parser)
// Complex/nested param defaults and all statements live inside workflow TrackTx.
// ============================================================================

params.reports_plots    = params.reports_plots    ?: 0
params.conda_norm       = params.conda_norm       ?: null
params.conda_pol        = params.conda_pol        ?: null

// ============================================================================
// TOP-LEVEL HELPER FUNCTION
// Must be declared at top level (not inside a workflow) in strict syntax.
// ============================================================================

// Resolve local FASTQ path: use as-is if absolute, else relative to projectDir
def resolveLocalPath(path) {
  if (!path?.trim()) return null
  def p = path.trim()
  def f = new File(p)
  return f.isAbsolute() ? file(p) : file("${projectDir}/${p}")
}

// ============================================================================
// MODULE IMPORTS
// Paths inlined (def MOD = ... was a top-level statement, not allowed in strict)
// ============================================================================

include { download_genome_annotations                               } from "${projectDir}/modules/01_download_genome_annotations.nf"
include { download_sra_samples                                      } from "${projectDir}/modules/02_download_sra_samples.nf"
include { preprocess_and_quality_filter_reads                       } from "${projectDir}/modules/03_preprocess_and_quality_filter_reads.nf"
include { download_genome_and_build_alignment_index as build_index  } from "${projectDir}/modules/04_download_genome_and_build_alignment_index.nf"
include { download_genome_and_build_alignment_index as spike_index  } from "${projectDir}/modules/04_download_genome_and_build_alignment_index.nf"
include { align_reads_to_genome                                     } from "${projectDir}/modules/05_align_reads_to_genome.nf"
include { check_and_merge_replicates                                } from "${projectDir}/modules/05b_check_and_merge_replicates.nf"
include { generate_coverage_tracks                                  } from "${projectDir}/modules/06_generate_coverage_tracks.nf"
include { quantify_reads_per_gene                                   } from "${projectDir}/modules/07_quantify_reads_per_gene.nf"
include { normalize_coverage_tracks                                 } from "${projectDir}/modules/08_normalize_coverage_tracks.nf"
include { detect_divergent_transcription                            } from "${projectDir}/modules/09_detect_divergent_transcription.nf"
include { assign_signal_to_functional_regions                       } from "${projectDir}/modules/10_assign_signal_to_functional_regions.nf"
include { calculate_polymerase_occupancy_metrics                    } from "${projectDir}/modules/11_calculate_polymerase_occupancy_metrics.nf"
include { summarize_polymerase_metrics                              } from "${projectDir}/modules/12_summarize_polymerase_metrics.nf"
include { quality_control_aligned_reads                             } from "${projectDir}/modules/13_quality_control_aligned_reads.nf"
include { generate_per_sample_reports                               } from "${projectDir}/modules/14_generate_per_sample_reports.nf"
include { combine_reports_into_cohort                               } from "${projectDir}/modules/15_combine_reports_into_cohort.nf"
include { cohort_qc_and_viz                                         } from "${projectDir}/modules/16_cohort_qc_and_viz.nf"

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {

  // ── Nested param defaults (complex/conditional — must live inside workflow) ─

  def qcParams  = params.qc  instanceof Map ? params.qc  : [:]
  def polParams = params.pol instanceof Map ? params.pol : [:]

  // ── Help message ───────────────────────────────────────────────────────────

  if (params.help) {
    log.info """
    ════════════════════════════════════════════════════════════════════════
    TrackTx Pipeline — PRO-seq Nascent RNA Analysis
    ════════════════════════════════════════════════════════════════════════

    Usage:
      nextflow run main.nf -profile docker \\
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

    For detailed documentation, see README.md
    ════════════════════════════════════════════════════════════════════════
    """.stripIndent()
    exit 0
  }

  // ── Startup banner ─────────────────────────────────────────────────────────

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

  // ── Parameter validation ───────────────────────────────────────────────────

  if (!params.output_dir) {
    error "PIPELINE | ERROR | Missing required parameter: --output_dir"
  }
  if (!params.samplesheet) {
    error "PIPELINE | ERROR | Missing required parameter: --samplesheet"
  }
  if (!params.reference_genome) {
    error "PIPELINE | ERROR | Missing required parameter: --reference_genome"
  }
  if (params.reference_genome == 'other' && !params.genome_fasta) {
    error "PIPELINE | ERROR | When reference_genome=other, must provide --genome_fasta"
  }
  if (params.spikein_genome == 'other' && !params.spikein_fasta) {
    error "PIPELINE | ERROR | When spikein_genome=other, must provide --spikein_fasta"
  }

  // Validate samplesheet exists (resolve relative paths from projectDir)
  def samplesheetPath = params.samplesheet?.trim()
  def samplesheetFile = (samplesheetPath && new File(samplesheetPath).isAbsolute())
    ? file(samplesheetPath)
    : file("${projectDir}/${samplesheetPath ?: 'samplesheet.csv'}")
  if (!samplesheetFile.exists()) {
    error "PIPELINE | ERROR | Samplesheet not found: ${params.samplesheet}"
  }

  def samplesheetLines = samplesheetFile.readLines()
  def dataRows         = samplesheetLines.drop(1).findAll { it.trim() }
  if (dataRows.isEmpty()) {
    error "PIPELINE | ERROR | Samplesheet has no data rows. Expected format: sample,condition,timepoint,replicate,file1,file2"
  }

  // Fail-fast: validate input files exist for local samples
  if (params.sample_source != 'srr') {
    def header   = samplesheetLines[0]?.split(',')?.collect { it.trim() }
    def sampleIdx = header?.findIndexOf { it?.toLowerCase() == 'sample' }
    def file1Idx  = header?.findIndexOf { it?.toLowerCase() == 'file1' }
    def file2Idx  = header?.findIndexOf { it?.toLowerCase() == 'file2' }
    if (sampleIdx != null && sampleIdx >= 0 && file1Idx != null && file1Idx >= 0) {
      def missing = []
      dataRows.eachWithIndex { line, i ->
        def cols   = line.split(',', -1).collect { it?.trim() }
        def sample = cols.size() > sampleIdx ? cols[sampleIdx] : ''
        def resolve = { p ->
          if (!p?.trim()) return null
          def s = p.trim()
          def f = new File(s)
          return f.isAbsolute() ? f : new File(projectDir.toString(), s)
        }
        def f1 = cols.size() > file1Idx ? resolve(cols[file1Idx]) : null
        def f2 = (params.paired_end && file2Idx != null && file2Idx >= 0 && cols.size() > file2Idx)
          ? resolve(cols[file2Idx]) : null
        if (f1 && !f1.exists()) missing << "${sample}: file1 not found: ${f1}"
        if (f2 && !f2.exists()) missing << "${sample}: file2 not found: ${f2}"
      }
      if (!missing.isEmpty()) {
        error """PIPELINE | ERROR | Input files not found (${missing.size()} issue(s)).
Fix paths in samplesheet or ensure files exist before running.

Missing/not found:
${missing.take(10).join('\n')}${missing.size() > 10 ? '\n... and ' + (missing.size() - 10) + ' more' : ''}

Expected format: sample,condition,timepoint,replicate,file1,file2
Paths are relative to: ${projectDir}"""
      }
    }
  }

  if (params.verbose) log.info "PIPELINE | VALIDATE | Parameter validation complete"
  if (params.verbose) log.info "PIPELINE | IMPORT | All modules loaded successfully"

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 1: Download Annotations
  // ══════════════════════════════════════════════════════════════════════════

  def assetsDir0  = params.assets_dir ?: "${projectDir}/assets"
  def noGtfPath0  = "${assetsDir0}/NO_GTF"
  if (!new File(noGtfPath0).exists()) new File(noGtfPath0).text = ''
  def customAnnotationFile = (params.gtf_path?.trim())
    ? file(params.gtf_path, checkIfExists: true)
    : file(noGtfPath0)

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 1 | Download Genome Annotations"
    log.info "-".multiply(80)
  }

  download_genome_annotations(channel.value(customAnnotationFile))

  def gtf_ch   = download_genome_annotations.out.gtf
  def genes_ch = download_genome_annotations.out.genes
  def tss_ch   = download_genome_annotations.out.tss
  def tes_ch   = download_genome_annotations.out.tes

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 2: Parse Samplesheet
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 2 | Parse Samplesheet"
    log.info "-".multiply(80)
    log.info "STEP 2 | INPUT | File: ${params.samplesheet}"
  }

  def samples_ch = channel
    .fromPath(samplesheetFile)
    .splitCsv(header: true)
    .map { row ->
      if (!row.sample || !row.file1) {
        error "STEP 2 | ERROR | Samplesheet row missing 'sample' or 'file1': ${row}"
      }
      if (params.sample_source != 'srr' && params.paired_end && !row.file2) {
        error "STEP 2 | ERROR | Paired-end mode requires 'file2' for sample ${row.sample}"
      }
      def sample_id = "${row.sample}_r${row.replicate ?: 1}"
      def reads = (params.sample_source == 'srr')
        ? [row.file1.trim(), row.file2?.trim()]
        : (params.paired_end
            ? [resolveLocalPath(row.file1), resolveLocalPath(row.file2)]
            : [resolveLocalPath(row.file1)])
      tuple(
        sample_id,
        reads,
        (row.treatment ?: row.condition ?: row.sample) ?: '',
        row.timepoint ?: '',
        row.replicate ?: 1
      )
    }
    .tap { samples_for_count }

  samples_for_count.count().subscribe { count ->
    if (count == 0) {
      error "PIPELINE | ERROR | No samples parsed from samplesheet. Check column names (sample, file1), delimiter (comma), and encoding (UTF-8)."
    }
    if (params.verbose) log.info "STEP 2 | COMPLETE | Loaded ${count} samples from samplesheet"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 3: Download SRR Data (optional)
  // ══════════════════════════════════════════════════════════════════════════

  def prepared_input_ch    = null
  def preexisting_clean_ch = channel.empty()

  if (params.sample_source == 'srr') {
    if (params.verbose) {
      log.info "-".multiply(80)
      log.info "STEP 3 | Download SRR Data from NCBI"
      log.info "-".multiply(80)
    }

    def isPairedEnd = params.paired_end ? true : false
    def outputDir   = params.output_dir.toString()

    if (params.verbose) log.info "STEP 3 | CHECK | Scanning for pre-existing trimmed FASTQs in ${outputDir}/01_trimmed_fastq/"

    def samples_branched = samples_ch.branch { item ->
      def sid        = item[0]
      def r1_trimmed = file("${outputDir}/01_trimmed_fastq/${sid}/final_R1.fastq")
      trimmed_exists: r1_trimmed.exists() && r1_trimmed.size() > 0
      needs_download: true
    }

    preexisting_clean_ch = samples_branched.trimmed_exists.map { sid, reads, c, t, r ->
      def r1 = file("${outputDir}/01_trimmed_fastq/${sid}/final_R1.fastq")
      def r2 = file("${outputDir}/01_trimmed_fastq/${sid}/final_R2.fastq")
      if (params.verbose) log.info "STEP 3 | SKIP | ${sid}: trimmed FASTQs found in results — skipping download and preprocessing"
      tuple(sid, r1, r2, c, t, r)
    }

    download_sra_samples(
      samples_ch.map { sid, reads, c, t, r ->
        tuple(sid, reads[0], c, t, r, isPairedEnd)
      }
    )
    
    prepared_input_ch = download_sra_samples.out[0].map { sid, fq1, fq2, c, t, r ->
      tuple(sid, [file(fq1), file(fq2)].findAll(), c, t, r)
    }

    prepared_input_ch.subscribe { sid, _reads, _cond, _time, _rep ->
      if (params.verbose) log.info "STEP 3 | DOWNLOAD | ${sid} complete"
    }
    

  } else {
    if (params.verbose) {
      log.info "-".multiply(80)
      log.info "STEP 3 | Using Local FASTQ Files"
      log.info "-".multiply(80)
    }
    prepared_input_ch = samples_ch
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 4: Preprocess and Clean Input Reads
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 4 | Preprocess Input Reads"
    log.info "-".multiply(80)
    log.info "STEP 4 | CONFIG | Mode: ${params.paired_end ? 'Paired-end' : 'Single-end'}"
  }

  def (preprocessed_clean_ch, fastqc_ch) = preprocess_and_quality_filter_reads(
    prepared_input_ch,
    channel.value(params.paired_end ? 'PE' : 'SE')
  )

  def clean_fastq_ch = preprocessed_clean_ch.mix(preexisting_clean_ch)

  clean_fastq_ch.subscribe { sid, _r1, _r2opt, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 4 | CLEAN | ${sid} preprocessing complete"
  }

  // Create sentinel files for SE mode / missing tracks
  def assetsDir       = params.assets_dir ?: "${projectDir}/assets"
  new File(assetsDir).mkdirs()

  def noR2Path        = "${assetsDir}/NO_R2.fastq"
  def noBGPath        = "${assetsDir}/EMPTY.bedgraph"
  def noBGPosPath     = "${assetsDir}/EMPTY_POS.bedgraph"
  def noBGNegPath     = "${assetsDir}/EMPTY_NEG.bedgraph"
  def noBG5pPosPath   = "${assetsDir}/EMPTY_5P_POS.bedgraph"
  def noBG5pNegPath   = "${assetsDir}/EMPTY_5P_NEG.bedgraph"
  def noBGAm5pPosPath = "${assetsDir}/EMPTY_AM5P_POS.bedgraph"
  def noBGAm5pNegPath = "${assetsDir}/EMPTY_AM5P_NEG.bedgraph"
  def noSpikeFaPath   = "${assetsDir}/EMPTY_SPIKE.fa"
  def noSpikeIdxPath  = "${assetsDir}/EMPTY_SPIKE_INDEX.fa"

  [noR2Path, noBGPath, noBGPosPath, noBGNegPath,
   noBG5pPosPath, noBG5pNegPath, noBGAm5pPosPath, noBGAm5pNegPath,
   noSpikeFaPath, noSpikeIdxPath].each { path ->
    if (!new File(path).exists()) new File(path).text = ''
  }
  new File(noSpikeFaPath).text  = ">none\nN\n"
  new File(noSpikeIdxPath).text = ">none_index\nN\n"

  def clean_fastq_with_r2 = clean_fastq_ch.map { sid, r1, r2opt, c, t, r ->
    def r2_file = (r2opt && file(r2opt).exists() && file(r2opt).size() > 0)
      ? file(r2opt)
      : file(noR2Path)
    tuple(sid, r1, r2_file, c, t, r)
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 5: Build Genome Indices
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 5 | Build Bowtie2 Indices"
    log.info "-".multiply(80)
  }

  def reference_fa = params.reference_genome == 'other'
    ? file(params.genome_fasta)
    : file("${projectDir}/genomes/${params.reference_genome}.fa")

  if (params.verbose) log.info "STEP 5 | INDEX | Primary genome: ${params.reference_genome}"

  build_index(
    channel.value(tuple(
      params.reference_genome == 'other' ? 'custom' : params.reference_genome,
      params.reference_genome == 'other' ? 'custom' : 'ucsc',
      reference_fa
    ))
  )

  def ref_meta_ch = build_index.out.ref_meta
  def ref_idx_ch  = build_index.out.index_files

  if (params.verbose) log.info "STEP 5 | INDEX | Primary genome index built"

  def noSpikeFa  = file(noSpikeFaPath)
  def noSpikeIdx = file(noSpikeIdxPath)
  def spike_meta_ch = channel.value(tuple('none', 'none', noSpikeFa))
  def spike_idx_ch  = channel.value(noSpikeIdx)

  if (params.spikein_genome && params.spikein_genome != 'None') {
    if (params.verbose) log.info "STEP 5 | INDEX | Spike-in genome: ${params.spikein_genome}"

    def spike_fa = params.spikein_genome == 'other'
      ? file(params.spikein_fasta)
      : file("${projectDir}/genomes/${params.spikein_genome}.fa")

    spike_index(
      channel.value(tuple(
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
    log.info "-".multiply(80)
    log.info "STEP 6 | Align Reads with Bowtie2"
    log.info "-".multiply(80)
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

  def aligned_ch = align_reads_to_genome.out[0]

  aligned_ch.subscribe { sid, _bam, _allbam, _spike, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 6 | ALIGN | ${sid} alignment complete"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 6b: Replicate Concordance Check & BAM Merging (optional)
  // ══════════════════════════════════════════════════════════════════════════

  if (params.replicates?.merge == true) {
    if (params.verbose) {
      log.info "-".multiply(80)
      log.info "STEP 6b | Replicate Concordance Check & BAM Merging"
      log.info "-".multiply(80)
      log.info "STEP 6b | CONFIG | Method: ${params.replicates?.concordance_method ?: 'pearson'}"
      log.info "STEP 6b | CONFIG | Threshold: ${params.replicates?.concordance_min ?: 0.9}"
    }

    def merge_input_ch = aligned_ch
      .map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
        tuple(c, t, sid, filt_bam, all_bam, spike_bam ?: file('NO_SPIKE'))
      }
      .groupTuple(by: [0, 1])

    check_and_merge_replicates(merge_input_ch)

    def merged_aligned_ch = check_and_merge_replicates.out.merged_bams
      .map { cond, tpt, filt_bam, all_bam, spike_bam ->
        def merged_sid = "${cond}_${tpt}_merged".toString().replaceAll(/[^a-zA-Z0-9_-]/, '_')
        def spike_f = (spike_bam.name != 'NO_SPIKE' && spike_bam.size() > 0) ? spike_bam : null
        tuple(merged_sid, filt_bam, all_bam, spike_f, cond, tpt, 0)
      }

    def passthrough_aligned_ch = check_and_merge_replicates.out.passthrough_bams
      .map { cond, tpt, filt_bams, all_bams, spike_bams ->
        [filt_bams].flatten().withIndex().collect { bam, i ->
          def allbam   = [all_bams].flatten()[i]
          def spikebam = [spike_bams].flatten()[i]
          def sid = bam.name.replaceAll(/\.bam$/, '')
          tuple(sid, bam, allbam, spikebam, cond, tpt, i + 1)
        }
      }
      .flatMap()

    aligned_ch = merged_aligned_ch.mix(passthrough_aligned_ch)

    check_and_merge_replicates.out.concordance_report
      .collectFile(
        name:     'concordance_report.tsv',
        storeDir: "${params.output_dir}/02_alignments/_merged",
        keepHeader: true,
        newLine:    false
      )

    aligned_ch.subscribe { sid, _bam, _allbam, _spike, _cond, _time, _rep ->
      if (params.verbose) log.info "STEP 6b | MERGE | ${sid} ready for track generation"
    }
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 7: Generate Coverage Tracks
  // ══════════════════════════════════════════════════════════════════════════
  
  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 7 | Generate Coverage Tracks"
    log.info "-".multiply(80)
  }

  def genome_fa_ch = ref_meta_ch.map { id, prov, fa -> fa }

  def tracks_input_ch = aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
    tuple(sid, filt_bam, spike_bam, c, t, r)
  }

  def allmap_bam_ch = aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
    all_bam
  }

  generate_coverage_tracks(
    tracks_input_ch,
    genome_fa_ch,
    channel.value(params.paired_end ? 'true' : 'false'),
    allmap_bam_ch
  )

  def bw3p_pair_ch     = generate_coverage_tracks.out.bw3p_pair
  def bw5p_pair_ch     = generate_coverage_tracks.out.bw5p_pair
  def allmap3p_pair_ch = generate_coverage_tracks.out.allmap3p_pair
  def allmap5p_pair_ch = generate_coverage_tracks.out.allmap5p_pair
  def tracks_ch        = generate_coverage_tracks.out.track_tuple
  def dedup_stats_ch   = generate_coverage_tracks.out.dedup_stats

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
    log.info "-".multiply(80)
    log.info "STEP 8 | Collect Read Counts"
    log.info "-".multiply(80)
  }

  def counts_tsvs = quantify_reads_per_gene(
    aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
      tuple(sid, filt_bam, all_bam, spike_bam ?: '-', c, t, r)
    }
  )

  def counts_master = counts_tsvs.counts
    .map { sid, tsv, c, t, r -> tsv }
    .collectFile(
      name:      'counts_master.tsv',
      storeDir:  "${params.output_dir}/04_counts",
      keepHeader: true,
      newLine:    true
    )

  counts_master.subscribe { f ->
    if (params.verbose) log.info "STEP 8 | COMPLETE | Master counts file: ${f}"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 9: Normalize Tracks
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 9 | Normalize Tracks"
    log.info "-".multiply(80)
  }

  def norm_main_kv = tracks_ch.map { sid, bam, spk, p3, n3, p5, n5, c, t, r ->
    def p3_file = file(p3)
    def n3_file = file(n3)
    def p5_file = (p5 && file(p5).exists() && file(p5).size() > 0) ? file(p5) : file(noBG5pPosPath)
    def n5_file = (n5 && file(n5).exists() && file(n5).size() > 0) ? file(n5) : file(noBG5pNegPath)
    tuple(sid, tuple(p3_file, n3_file, p5_file, n5_file, c, t, r))
  }

  def norm_allmap3_kv = allmap3p_pair_ch.map { sid, ap3, an3, bwp, bwn, c, t, r ->
    tuple(sid, tuple(file(ap3), file(an3)))
  }

  def norm_allmap5_kv = allmap5p_pair_ch.map { sid, ap5, an5, bwp, bwn, c, t, r ->
    def ap5_file = (ap5 && file(ap5).exists() && file(ap5).size() > 0) ? file(ap5) : file(noBGAm5pPosPath)
    def an5_file = (an5 && file(an5).exists() && file(an5).size() > 0) ? file(an5) : file(noBGAm5pNegPath)
    tuple(sid, tuple(ap5_file, an5_file))
  }

  def norm_input_ch = norm_main_kv
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

  normalize_coverage_tracks(norm_input_ch, genome_fa_ch, tes_ch)

  def norm_tracks_ch  = normalize_coverage_tracks.out.norm_tuple
  def norm_factors_ch = norm_tracks_ch.map { sid, p3, n3, nf, c, t, r ->
    tuple(sid, nf, c, t, r)
  }

  norm_tracks_ch.subscribe { sid, _p3, _n3, _nf, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 9 | NORMALIZE | ${sid} normalization complete"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 10: Detect Divergent Transcription
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "=".multiply(80)
    log.info "STEP 10 | Detect Divergent Transcription"
    log.info "=".multiply(80)
  }

  def divergent_input_ch = bw3p_pair_ch
    .map { sid, pos3_bg, neg3_bg, bwp, bwn, c, t, r ->
      tuple(sid, pos3_bg, neg3_bg, c, t, r)
    }

  divergent_input_ch.subscribe { sid, _pos, _neg, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 10 | INPUT | Sample ready for detection: ${sid}"
  }

  detect_divergent_transcription(
    divergent_input_ch,
    channel.value(params.advanced?.divergent_threshold ?: 'auto'),
    channel.value(params.advanced?.divergent_sum_thr ?: 'auto'),
    channel.value(params.advanced?.divergent_fdr ?: 0.08),
    channel.value(params.advanced?.divergent_nt_window ?: 1000),
    channel.value(params.advanced?.divergent_balance ?: 0.0),
    channel.value(params.advanced?.divergent_bin_gap ?: 100),
    channel.value(params.advanced?.divergent_calibration_percentile ?: 65.0),
    channel.value(params.advanced?.divergent_calibration_sum_multiplier ?: 1.5),
    channel.value(params.advanced?.divergent_calibration_background_lower ?: false),
    channel.value(params.advanced?.divergent_merge_gap ?: 150)
  )

  def divergent_tx_ch = detect_divergent_transcription.out.bed

  divergent_tx_ch.subscribe { sid, bed, _cond, _time, _rep ->
    def count = 0
    if (file(bed).exists() && file(bed).size() > 0) {
      count = file(bed).readLines().findAll { ln -> !ln.startsWith('#') }.size()
    }
    if (params.verbose) log.info "STEP 10 | COMPLETE | ${sid} → ${count} divergent regions detected"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 11: Score Divergent Sites (Enhancer vs Gene)
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 11 | Call Functional Regions"
    log.info "-".multiply(80)
  }

  def func_input_ch = divergent_tx_ch
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
      tuple(sid, div_bed, pos3_raw, neg3_raw, file(noBGPosPath), file(noBGNegPath), c, t, r)
    }

  assign_signal_to_functional_regions(
    func_input_ch,
    gtf_ch,
    channel.value(file("${projectDir}/bin/functional_regions.py")),
    genes_ch,
    tss_ch,
    tes_ch
  )

  def functional_regions_ch = assign_signal_to_functional_regions.out.main

  def functional_regions_bed_ch = functional_regions_ch.map { sid, bed, fsum, c, t, r ->
    tuple(sid, bed)
  }

  def functional_regions_sum_ch = functional_regions_ch.map { sid, bed, fsum, c, t, r ->
    tuple(sid, fsum, c, t, r)
  }

  functional_regions_ch.subscribe { sid, bed, _fsum, _cond, _time, _rep ->
    def count = 0
    if (file(bed).exists() && file(bed).size() > 0) {
      count = file(bed).readLines().findAll { ln -> !ln.startsWith('#') }.size()
    }
    if (params.verbose) log.info "STEP 12 | COMPLETE | ${sid} → ${count} functional regions annotated"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 12: Calculate Pol-II Metrics
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 12 | Calculate Pol-II Metrics"
    log.info "-".multiply(80)
  }

  def pol_input_ch = generate_coverage_tracks.out.bam_for_tracks
    .map { sid, bam, c, t, r ->
      tuple(sid, tuple(bam, c, t, r))
    }
    .join(functional_regions_bed_ch)
    .join(
      norm_tracks_ch.map { sid, pos3_cpm, neg3_cpm, factors, c, t, r ->
        def normDir     = "${params.output_dir}/05_normalized_tracks/${sid}"
        def pos3_sicpm  = file("${normDir}/sicpm/3p/${sid}.3p.pos.sicpm.bedgraph")
        def neg3_sicpm  = file("${normDir}/sicpm/3p/${sid}.3p.neg.sicpm.bedgraph")
        tuple(sid, tuple(pos3_cpm, neg3_cpm, pos3_sicpm, neg3_sicpm))
      }
    )
    .map { sid, bam_data, bed, norm_data ->
      def (bam, c, t, r) = bam_data
      def (pos_cpm, neg_cpm, pos_si, neg_si) = norm_data
      tuple(sid, bam, bed, pos_cpm, neg_cpm, pos_si, neg_si, c, t, r)
    }

  calculate_polymerase_occupancy_metrics(pol_input_ch, gtf_ch)

  def pol_gene_ch    = calculate_polymerase_occupancy_metrics.out.genes
  def pol_density_ch = calculate_polymerase_occupancy_metrics.out.density
  def pol_pausing_ch = calculate_polymerase_occupancy_metrics.out.pausing

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 13: Summarize Pol-II Metrics (Cohort-Level)
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 13 | Summarize Pol-II Metrics"
    log.info "-".multiply(80)
  }

  def pol_sorted = pol_gene_ch
    .toSortedList { a, b -> a[0] <=> b[0] }

  def samples_lines = pol_sorted
    .flatMap { sorted_list ->
      sorted_list.withIndex().collect { item, idx ->
        def (sid, genes, c, t, r) = item
        "${sid}\t${c ?: 'NA'}\t${t ?: 'NA'}\t${r ?: '1'}\tmetric_${idx + 1}"
      }
    }

  def samples_tsv = channel
    .of('sample_id\tcondition\ttimepoint\treplicate\tfile')
    .concat(samples_lines)
    .collectFile(
      name:    'samples.tsv',
      newLine: true
    )
    .map { tsv_file ->
      def clean_text = tsv_file.text
        .replace('﻿', '')
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
        def rawCount  = lines.size() - 1
        def rejectMsg = rejected ? "\n  Rejected rows (first 5): ${rejected.take(5).join('\n  ')}" : ''
        error """STEP 12 | ERROR | No valid samples in Pol-II aggregate TSV
Raw data rows: ${rawCount} | Valid after parse: ${clean_lines.size() - 1}${rejectMsg}
This usually means no samples reached calculate_polymerase_occupancy_metrics."""
      }

      tsv_file.text = clean_lines.join('\n') + '\n'
    if (params.verbose) log.info "STEP 14 | INPUT | Samples TSV prepared: ${clean_lines.size() - 1} samples"
      tsv_file
    }

  def pol_files_ch = pol_sorted
    .flatMap { sorted_list ->
      sorted_list.collect { sid, genes, c, t, r -> file(genes) }
    }

  summarize_polymerase_metrics(samples_tsv, pol_files_ch.collect())

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 14: Quality Control
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 14 | Quality Control Analysis"
    log.info "-".multiply(80)
  }

  def qc_input_ch = aligned_ch
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

  def qc_json_meta_ch = quality_control_aligned_reads.out.json_meta

  qc_json_meta_ch.subscribe { sid, _json, _cond, _time, _rep ->
    if (params.verbose) log.info "STEP 14 | COMPLETE | ${sid} QC analysis finished"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 15: Generate Per-Sample Reports
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 15 | Generate Per-Sample Reports"
    log.info "-".multiply(80)
  }

  def resolvePath = { path ->
    def f = file(path)
    f.exists() ? f.toString() : ''
  }

  def report_input_ch = aligned_ch
    .map { sid, bam, all, spike, c, t, r ->
      tuple(sid, tuple(c ?: 'NA', t ?: 'NA', r ?: '1'))
    }
    .join(divergent_tx_ch.map  { sid, bed, c, t, r   -> tuple(sid, bed) })
    .join(functional_regions_sum_ch.map { sid, fsum, c, t, r -> tuple(sid, fsum) })
    .join(pol_density_ch.map   { sid, dens, c, t, r  -> tuple(sid, dens) })
    .join(pol_pausing_ch.map   { sid, paus, c, t, r  -> tuple(sid, paus) })
    .join(norm_factors_ch.map  { sid, nf, c, t, r    -> tuple(sid, nf) })
    .join(dedup_stats_ch.map   { sid, dedup, c, t, r -> tuple(sid, dedup) })
    .join(qc_json_meta_ch.map  { sid, qc, c, t, r    -> tuple(sid, qc) })
    .map { sid, meta, div_bed, fsum, dens, paus, norm, dedup, qc ->
      def (c, t, r) = meta
      def normDir   = "${params.output_dir}/05_normalized_tracks/${sid}"
      def tracksDir = "${params.output_dir}/03_genome_tracks/${sid}"
      def bw_pos3          = resolvePath("${normDir}/cpm/3p/${sid}.3p.pos.cpm.bw")
      def bw_neg3          = resolvePath("${normDir}/cpm/3p/${sid}.3p.neg.cpm.bw")
      def bw_allmap_pos3   = resolvePath("${normDir}/cpm/3p/${sid}.allMap.3p.pos.cpm.bw")
      def bw_allmap_neg3   = resolvePath("${normDir}/cpm/3p/${sid}.allMap.3p.neg.cpm.bw")
      def raw_allmap_pos3  = resolvePath("${tracksDir}/3p/${sid}.allMap.3p.pos.bedgraph")
      def raw_allmap_neg3  = resolvePath("${tracksDir}/3p/${sid}.allMap.3p.neg.bedgraph")
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

  def report_json_ch = generate_per_sample_reports.out[2]

  report_json_ch.subscribe { json_file ->
    if (params.verbose) log.info "STEP 15 | COMPLETE | Report generated: ${json_file.name}"
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 16: Cohort QC — MultiQC, deepTools PCA/Correlation, IGV, Run-on
  // (runs first so its outputs feed into the final landing page in Step 17)
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 16 | Cohort QC — MultiQC / deepTools / IGV / Run-on Efficiency"
    log.info "-".multiply(80)
  }

  // Collect per-sample BigWig and bedGraph paths in a consistent order.
  // We join all track channels on sample_id to guarantee ordering alignment,
  // then collect() into lists for the single cohort-level process invocation.

  // bw3p_pair_ch tuple: (sid, p3bg, n3bg, bwp, bwn, c, t, r)
  // After joins the cohort_tracks_ch tuple is:
  //   (sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan)
  def cohort_tracks_ch = bw3p_pair_ch
    .join(bw5p_pair_ch.map     { sid, p5bg, n5bg, bwp5, bwn5, c, t, r -> tuple(sid, p5bg, n5bg) })
    .join(allmap3p_pair_ch.map { sid, ap3,  an3,  bwap, bwan, c, t, r -> tuple(sid, bwap, bwan) })

  // Gather across all samples into sorted lists (toSortedList by sample_id)
  def cohort_bw_pos3    = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> bwp.toString()   }.toSortedList()
  def cohort_bw_neg3    = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> bwn.toString()   }.toSortedList()
  def cohort_bw_ampos3  = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> bwap.toString()  }.toSortedList()
  def cohort_bw_amneg3  = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> bwan.toString()  }.toSortedList()
  def cohort_sample_ids = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> sid              }.toSortedList()
  def cohort_conditions = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> c ?: 'unknown'  }.toSortedList()
  def cohort_pos3_bg    = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> p3bg.toString()  }.toSortedList()
  def cohort_neg3_bg    = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> n3bg.toString()  }.toSortedList()
  def cohort_pos5_bg    = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> p5bg.toString()  }.toSortedList()
  def cohort_neg5_bg    = cohort_tracks_ch.map { sid, p3bg, n3bg, bwp, bwn, c, t, r, p5bg, n5bg, bwap, bwan -> n5bg.toString()  }.toSortedList()

  // Collect all QC log files staged into a single directory for MultiQC.
  // Includes: bowtie2 logs, flagstats, trimming logs — all already emitted
  // from earlier modules as path outputs that land in the work directory.
  def multiqc_logs_ch = channel.empty()
    .mix(align_reads_to_genome.out.align_logs.flatten())
    .mix(align_reads_to_genome.out.flagstats.flatten())
    .collect()

  cohort_qc_and_viz(
    multiqc_logs_ch,
    cohort_bw_pos3,
    cohort_bw_neg3,
    cohort_bw_ampos3,
    cohort_bw_amneg3,
    cohort_sample_ids,
    cohort_conditions,
    cohort_pos3_bg,
    cohort_neg3_bg,
    cohort_pos5_bg,
    cohort_neg5_bg,
    genes_ch
  )

  if (params.verbose) {
    cohort_qc_and_viz.out.runon_efficiency.subscribe { tsv ->
      log.info "STEP 16 | COMPLETE | Run-on efficiency: ${tsv.name}"
    }
  }

  // ══════════════════════════════════════════════════════════════════════════
  // STEP 17: Combine Reports + Landing Page  (final step — runs last)
  // Module 16 outputs feed into the landing page so it is the definitive
  // entry point for the full run.
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    log.info "-".multiply(80)
    log.info "STEP 17 | Combine Reports + Landing Page"
    log.info "-".multiply(80)
  }

  def per_sample_reports = report_json_ch
    .toSortedList { a, b -> a.name <=> b.name }

  def concordance_for_cohort = (params.replicates?.merge == true)
    ? channel.fromPath("${params.output_dir}/02_alignments/_merged/concordance_report.tsv", checkIfExists: false)
             .ifEmpty(file("NO_CONCORDANCE"))
    : channel.value(file("NO_CONCORDANCE"))

  // Optional outputs from module 16 — use sentinel file when not produced
  // (e.g. deepTools plots require ≥2 samples; MultiQC may not be installed)
  def qc_multiqc_html  = cohort_qc_and_viz.out.multiqc_html .ifEmpty(file("NO_FILE"))
  def qc_igv_session   = cohort_qc_and_viz.out.igv_session
  def qc_runon_tsv     = cohort_qc_and_viz.out.runon_efficiency
  def qc_pca_plot      = cohort_qc_and_viz.out.pca_plot     .ifEmpty(file("NO_FILE"))
  def qc_corr_heatmap  = cohort_qc_and_viz.out.corr_heatmap .ifEmpty(file("NO_FILE"))

  combine_reports_into_cohort(
    per_sample_reports,
    concordance_for_cohort,
    qc_multiqc_html,
    qc_igv_session,
    qc_runon_tsv,
    qc_pca_plot,
    qc_corr_heatmap
  )

  // ══════════════════════════════════════════════════════════════════════════
  // PIPELINE COMPLETE
  // ══════════════════════════════════════════════════════════════════════════

  if (params.verbose) {
    def TIMESTAMP_END = new Date().format("yyyy-MM-dd HH:mm:ss")
    log.info "=".multiply(80)
    log.info "PIPELINE COMPLETE | TrackTx Analysis Finished Successfully"
    log.info "=".multiply(80)
    log.info "Started:   ${TIMESTAMP_START}"
    log.info "Finished:  ${TIMESTAMP_END}"
    log.info "Results:   ${params.output_dir}"
    log.info "=".multiply(80)
    log.info ""
    log.info "Next Steps:"
    log.info "  • MultiQC report:        ${params.output_dir}/12_cohort_qc/multiqc/multiqc_report.html"
    log.info "  • deepTools PCA:         ${params.output_dir}/12_cohort_qc/deeptools/pca_plot.pdf"
    log.info "  • Correlation heatmap:   ${params.output_dir}/12_cohort_qc/deeptools/correlation_heatmap.pdf"
    log.info "  • IGV session:           ${params.output_dir}/12_cohort_qc/igv_session.xml"
    log.info "  • Run-on efficiency:     ${params.output_dir}/12_cohort_qc/runon_efficiency.tsv"
    log.info "  • Per-sample QC:         ${params.output_dir}/10_qc/"
    log.info "  • HTML reports:          ${params.output_dir}/11_reports/"
    log.info "  • Divergent tx:          ${params.output_dir}/06_divergent_tx/"
    log.info "  • Pol-II metrics:        ${params.output_dir}/08_pol_metrics/"
    log.info "=".multiply(80)
  }
}
