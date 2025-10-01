// ╔═══════════════════════════════════════════════════════════════════════╗
/* ║ TrackTx — Nascent RNA processing & functional-region discovery (DSL2) ║
   ║                                                                       ║
   ║ Inputs  : FASTQ files (samplesheet) or SRR accessions                 ║
   ║ Outputs : CPM/siCPM bigWigs (no RPKM), divergent-Tx BED (allMap ±),   ║
   ║           functional regions (siCPM→CPM), Pol-II metrics, QC, reports ║ */
// ╚═══════════════════════════════════════════════════════════════════════╝

nextflow.enable.dsl = 2

// ───────────────────────────── 0) Parameter Defaults ───────────────────────────
// Additional parameter defaults (core ones are in nextflow.config)

params.reports_plots     = (params.reports_plots == null) ? 0 : params.reports_plots
params.conda_norm        = params.conda_norm ?: null
params.conda_pol2        = params.conda_pol2 ?: null

// QC defaults (used by qc_pol2_tracktx)
params.qc                = (params.qc instanceof Map) ? params.qc : [:]
params.qc.tool           = (params.qc.tool ?: 'fastp')
params.qc.mapq           = (params.qc.mapq          ?: 10)
params.qc.dedup          = (params.qc.dedup == null) ? true : params.qc.dedup
params.qc.depth_max_cov  = (params.qc.depth_max_cov ?: 0)

// Pol2 aggregate defaults
params.pol2              = (params.pol2 instanceof Map) ? params.pol2 : [:]
params.pol2.top_n        = (params.pol2.top_n ?: 100)
params.pol2.plots        = (params.pol2.plots == null) ? true : params.pol2.plots

// ───────────────────────────── 1) Schema & Help ───────────────────────────
// Schema validation disabled due to filesystem symlink issues with nf-schema plugin
// include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Simple help handling without nf-schema plugin
if( params.help ){
  log.info """
TrackTx Pipeline - Nascent RNA processing & functional-region discovery

Usage:
  nextflow run main.nf -entry TrackTx -profile docker -with-report -with-timeline -with-dag

Required parameters:
  --samplesheet <CSV>           Path to samplesheet CSV file
  --reference_genome <GENOME>   Reference genome (hg38, mm39, etc.)
  
Optional parameters:
  --output_dir <DIR>            Output directory (default: ./results)
  --paired_end                  Enable paired-end processing
  --debug                       Enable debug mode
  
For more details, see the README.md file.
"""
  exit 0
}

// Log summarized params (simplified)
log.info "Launching TrackTx pipeline..."
log.info "Output directory: ${params.output_dir}"
log.info "Samplesheet: ${params.samplesheet}"
log.info "Reference genome: ${params.reference_genome}"

// ───────────────────────── 0b) Early parameter validation ────────────────
if( !params.output_dir ) error "Missing: --output_dir"
if( !params.samplesheet ) error "Missing: --samplesheet"
if( !params.reference_genome ) error "Missing: --reference_genome"
if( params.reference_genome == 'other' && !params.genome_fasta )
  error "When reference_genome=other, provide --genome_fasta"
if( params.spikein_genome == 'other' && !params.spikein_fasta )
  error "When spikein_genome=other, provide --spikein_fasta"

// ───────────────────────── 1) Channel Shapes (cheat sheet) ────────────────
// (unchanged comments)

// ─────────────────────────── 2) Module Imports ────────────────────────────
def MOD = "${projectDir}/modules"

include { download_gtf                         } from "${MOD}/01_download_gtf.nf"
include { download_srr                         } from "${MOD}/02_download_srr.nf"
include { prepare_input                        } from "${MOD}/03_prepare_input.nf"
include { fetch_and_build_index as build_index } from "${MOD}/04_fetch_and_build_index.nf"
include { fetch_and_build_index as spike_index } from "${MOD}/04_fetch_and_build_index.nf"
include { run_alignment                        } from "${MOD}/05_run_alignment.nf"
include { generate_tracks                      } from "${MOD}/06_generate_tracks.nf"
include { collect_counts                       } from "${MOD}/07_collect_counts.nf"
include { normalize_tracks                     } from "${MOD}/08_normalize_tracks.nf"
include { detect_divergent_tx                  } from "${MOD}/09_detect_divergent_tx.nf"
include { call_functional_regions              } from "${MOD}/10_call_functional_regions.nf"
include { calculate_pol2_metrics               } from "${MOD}/11_calculate_pol2_metrics.nf"
include { summarize_pol2_metrics               } from "${MOD}/12_summarize_pol2_metrics.nf"
include { qc_pol2_tracktx                      } from "${MOD}/13_qc_pol2_tracktx.nf"
include { generate_reports                     } from "${MOD}/14_generate_reports.nf"
include { combine_reports                      } from "${MOD}/15_combine_reports.nf"

// ───────────────────────────────── 3) Workflow ────────────────────────────
workflow TrackTx {

  // ───────── 3.1 Annotation ─────────
  download_gtf()
  def gtf_ch   = download_gtf.out[0]
  def genes_ch = download_gtf.out.genes
  def tss_ch   = download_gtf.out.tss
  def tes_ch   = download_gtf.out.tes

  // ───────── 3.2 Samplesheet ────────
  if( !params.samplesheet ) error "Missing: --samplesheet <CSV>"

  def samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      if (!row.sample || !row.file1)
        error "Samplesheet row missing 'sample' or 'file1': ${row}"
      if (params.paired_end && !row.file2)
        error "paired_end=true but 'file2' empty for sample ${row.sample}"

      def sid   = "${row.sample}_r${row.replicate ?: 1}"
      def reads =
        (params.sample_source == 'srr')
          ? [ row.file1.trim(), row.file2?.trim() ]
          : (params.paired_end ? [ file("${projectDir}/${row.file1}"), file("${projectDir}/${row.file2}") ] : [ file("${projectDir}/${row.file1}") ])
      // Map condition to the 'sample' column to align grouping and siCPM control selection
      tuple(sid, reads, row.sample, row.timepoint, row.replicate ?: 1)
    }

  // ───────── 3.3 SRR fetch (opt) ────
  def prepared_input_ch =
    (params.sample_source == 'srr')
      ? download_srr(
          samples_ch.map { sid, reads, c, t, r ->
            tuple(sid, reads[0], c, t, r, params.paired_end)
          }
        ).map { sid, fq1, fq2, c, t, r ->
          tuple(sid, [ file(fq1), file(fq2) ].findAll(), c, t, r)
        }
      : samples_ch

  // ───────── 3.4 Pre-process ────────
  def (clean_fastq_ch, fastqc_ch) =
    prepare_input(prepared_input_ch, Channel.value(params.paired_end ? 'PE' : 'SE'))

  // ───────── 3.4b Sentinel for SE (inputs cannot be optional) ────────────
  // Create a stable empty R2 path and inject it when R2 is absent.
  def assetsDir  = params.assets_dir ?: "${projectDir}/assets"
  new File(assetsDir).mkdirs()
  def noR2Path   = "${assetsDir}/NO_R2.fastq"
  def noBGPath   = "${assetsDir}/EMPTY.bedgraph"
  if( !new File(noR2Path).exists() ) new File(noR2Path).text = ''
  if( !new File(noBGPath).exists() ) new File(noBGPath).text = ''

  // cleaned_reads → (sid, r1, r2?, cond, tp, rep)  → ensure r2 path always present
  def clean_fastq_with_r2 =
    clean_fastq_ch.map { sid, r1, r2opt, c, t, r ->
      def r2 = (r2opt && file(r2opt).exists() && file(r2opt).size()>0) ? file(r2opt) : file(noR2Path)
      tuple(sid, r1, r2, c, t, r)
    }

  // ───────── 3.5 References ────────
  def reference_fa = params.reference_genome == 'other'
      ? file(params.genome_fasta)
      : file("${projectDir}/genomes/${params.reference_genome}.fa")

  def spike_fa = (params.spikein_genome in ['None', null]) ? null
               : (params.spikein_genome == 'other'
                    ? file(params.spikein_fasta)
                    : file("${projectDir}/genomes/${params.spikein_genome}.fa"))

  build_index(
    Channel.value(tuple(
      params.reference_genome == 'other' ? 'custom' : params.reference_genome,
      params.reference_genome == 'other' ? 'custom' : 'ucsc',
      reference_fa
    ))
  )
  def ref_meta_ch = build_index.out.ref_meta
  def ref_idx_ch  = build_index.out.index_files

  def spike_meta_ch
  def spike_idx_ch
  if (spike_fa) {
    spike_index(Channel.value(tuple(
      params.spikein_genome == 'other' ? 'custom_spike' : params.spikein_genome,
      params.spikein_genome == 'other' ? 'custom'       : 'ucsc',
      spike_fa
    )))
    spike_meta_ch = spike_index.out.ref_meta
    spike_idx_ch  = spike_index.out.index_files
  } else {
    spike_meta_ch = Channel.empty()
    spike_idx_ch  = Channel.empty()
  }

  // ───────── 3.6 Alignment ──────────
  // NOTE: run_alignment expects read2 as a required path; we pass NO_R2.fastq for SE.
  def (aligned_ch, align_logs_ch, flagstats_ch, idxstats_ch) =
    run_alignment(clean_fastq_with_r2, ref_meta_ch, ref_idx_ch, spike_meta_ch, spike_idx_ch)
  // aligned_ch: (sid, sample.bam, sample_allMap.bam, spikein.bam, c,t,r)

  // ───────── 3.7 Tracks ─────────────
  def genome_fa_ch = ref_meta_ch.map { _id, _prov, fa -> fa }

  def tracks_tuple_ch = aligned_ch.map { sid, filt_bam, _all_bam, spike_bam, c, t, r ->
    tuple(sid, filt_bam, spike_bam, c, t, r)
  }
  def allmap_bam_ch  = aligned_ch.map { sid, _filt_bam, all_bam, _spike_bam, _c, _t, _r -> all_bam }

  // output order: bw3p, bw5p, allmap3p, allmap5p, track_tuple, dedup_stats
  def (bw3p_pair_ch, bw5p_pair_ch, allmap3p_pair_ch, allmap5p_pair_ch, tracks_ch, dedup_stats_ch) =
    generate_tracks(
      tracks_tuple_ch,
      genome_fa_ch,
      Channel.value(params.paired_end ? 'true' : 'false'),
      allmap_bam_ch
    )

  if (params.debug) {
    tracks_ch.view { sid, fb, sb, p3, n3, p5, n5, c, t, r ->
      "DEBUG[tracks] ${sid} 3p(+/-)=${p3},${n3} 5p(+/-)=${p5 ?: '∅'},${n5 ?: '∅'}"
    }
    allmap3p_pair_ch.view { sid, ap, an, _bwp, _bwn, c,t,r -> "DEBUG[allMap3p] ${sid} allMap3p(+/-)=${ap},${an}" }
  }

  // ───────── 3.8 Read counts ────────
  def (counts_tsvs) = collect_counts(
    aligned_ch.map { sid, filt_bam, all_bam, spike_bam, c, t, r ->
      tuple(sid, filt_bam, all_bam, (spike_bam ?: '-'), c, t, r)
    }
  )

  def counts_master = counts_tsvs.collectFile(
    name: 'counts_master.tsv',
    storeDir: "${params.output_dir}/04_counts",
    keepHeader: true,
    newLine: true
  )

  // ───────── 3.9 Normalization ──────
  // Build full, fixed-arity tuples for normalize_tracks (no optional inputs).
  // Fill missing 5′ (main/allMap) with EMPTY.bedgraph placeholders.

  def norm_left_main =
    tracks_ch.map { sid, _bam, _spk, p3, n3, p5, n5, c, t, r ->
      def p5x = (p5 && file(p5).exists() && file(p5).size()>0) ? file(p5) : file(noBGPath)
      def n5x = (n5 && file(n5).exists() && file(n5).size()>0) ? file(n5) : file(noBGPath)
      tuple(sid, p3, n3, p5x, n5x, c, t, r)
    }

  def norm_left_allmap3 =
    allmap3p_pair_ch.map { sid, ap3, an3, _bwp, _bwn, c, t, r ->
      tuple(sid, ap3, an3, c, t, r)
    }

  def norm_left_allmap5 =
    aligned_ch.map { sid, _fb, _ab, _sb, c, t, r ->
      tuple(sid, file(noBGPath), file(noBGPath), c, t, r)
    }

  // Join into the exact argument list normalize_tracks expects
  def norm_input_ch =
    norm_left_main
      .join(norm_left_allmap3)
      .join(norm_left_allmap5)
      .combine(counts_master)   // singleton counts_master.tsv
      .combine(genes_ch)        // singleton genes (iface stability)
      .map { sid, p3, n3, p5, n5, c1, t1, r1,
                  ap3, an3, _c2, _t2, _r2,
                  ap5, an5, _c3, _t3, _r3,
                  cm, genes ->
        tuple(sid, p3, n3, p5, n5, ap3, an3, ap5, an5, c1, t1, r1, cm, genes)
      }

  normalize_tracks( norm_input_ch, genome_fa_ch )
  def norm_tracks_ch  = normalize_tracks.out.norm_tuple
  def norm_factors_ch = norm_tracks_ch.map { sid, _p3, _n3, nf, c, t, r -> tuple(sid, nf, c, t, r) }

  // ───────── 3.10 Divergent Tx (raw allMap 3′ ±) ─────────────────────────
  def divergent_input_ch =
    allmap3p_pair_ch.map { sid, ap3, an3, _bwp, _bwn, c, t, r ->
      tuple(sid, tuple(ap3, an3, c, t, r))
    }
    .groupTuple()
    .map { sid, values ->
      // values is a List of tuples: (ap3, an3, c, t, r)
      def pick = values.find { v ->
        def (ap, an, _c, _t, _r) = v
        (file(ap).exists() && file(ap).size()>0) || (file(an).exists() && file(an).size()>0)
      } ?: values[0]
      def (ap, an, c, t, r) = pick
      tuple(sid, ap, an, c, t, r)
    }

  detect_divergent_tx( divergent_input_ch )
  def divergent_tx_ch = detect_divergent_tx.out[0]

  // ───────── 3.11 Functional regions ─────────────────────────────────────
  def div_kv_fgr =
    divergent_tx_ch.map { sid, div_bed, c, t, r -> tuple(sid, tuple(div_bed, c, t, r)) }

  def tracks_3p_norm_kv_fgr =
    norm_tracks_ch.map { sid, pos3_cpm_bg, neg3_cpm_bg, _factors, c, t, r ->
      def pos3_sicpm_bg = file(pos3_cpm_bg.toString().replace('_pos3_cpm.bedgraph', '_pos3_sicpm.bedgraph'))
      def neg3_sicpm_bg = file(neg3_cpm_bg.toString().replace('_neg3_cpm.bedgraph', '_neg3_sicpm.bedgraph'))
      tuple(sid, tuple(pos3_cpm_bg, neg3_cpm_bg, pos3_sicpm_bg, neg3_sicpm_bg, c, t, r))
    }

  def func_input_ch =
    div_kv_fgr.join(tracks_3p_norm_kv_fgr)
      .map { sid, L, R ->
        def (div_bed, c1, t1, r1) = L
        def (pos3_cpm, neg3_cpm, pos3_sicpm, neg3_sicpm, _c2, _t2, _r2) = R
        tuple(sid, div_bed, pos3_cpm, neg3_cpm, pos3_sicpm, neg3_sicpm, c1, t1, r1)
      }

  call_functional_regions(
    func_input_ch,
    gtf_ch,
    Channel.value(file("${projectDir}/bin/functional_regions.py")),
    genes_ch, tss_ch, tes_ch
  )

  def functional_regions_ch     = call_functional_regions.out.main
  def functional_regions_bed_ch = functional_regions_ch.map { sid, bed,  _fsum, c, t, r -> tuple(sid, bed,  c, t, r) }
  def functional_regions_sum_ch = functional_regions_ch.map { sid, _bed, fsum,  c, t, r -> tuple(sid, fsum, c, t, r) }

  // ───────── 3.12 Pol II metrics (from normalized 3′) ────────────────────
  def pol2_left_kv =
    aligned_ch.map { sid, filt_bam, _all, _spike, c, t, r ->
      tuple(sid, [filt_bam, (c ?: ''), (t ?: ''), (r ?: 1)])
    }

  def bed_kv = functional_regions_bed_ch.map { sid, bed, _c, _t, _r -> tuple(sid, bed) }

  def pol2_left_with_bed_kv =
    pol2_left_kv.join(bed_kv)
      .map { it ->
        def (sid, L, bed) = it
        def (bam, cond, tp, rep) = L
        tuple(sid, [bam, bed, cond, tp, rep])
      }

  def pol2_norm_kv =
    norm_tracks_ch.map { sid, pos3_cpm_bg, neg3_cpm_bg, _factors, _c, _t, _r ->
      def pos3_si = file(pos3_cpm_bg.toString().replace('_pos3_cpm.bedgraph', '_pos3_sicpm.bedgraph'))
      def neg3_si = file(neg3_cpm_bg.toString().replace('_neg3_cpm.bedgraph', '_neg3_sicpm.bedgraph'))
      tuple(sid, [pos3_cpm_bg, neg3_cpm_bg, pos3_si, neg3_si])
    }

  def pol2_input_ch =
    pol2_left_with_bed_kv.join(pol2_norm_kv)
      .map { it ->
        def (sid, L, N) = it
        def (bam, bed, cond, tp, rep) = L
        def (pos_cpm, neg_cpm, pos_si, neg_si) = N
        tuple(sid, bam, bed, pos_cpm, neg_cpm, pos_si, neg_si, cond, tp, rep)
      }

  calculate_pol2_metrics(
    pol2_input_ch,
    gtf_ch
  )
  def pol2_gene_ch    = calculate_pol2_metrics.out.genes
  def pol2_density_ch = calculate_pol2_metrics.out.density
  def pol2_pausing_ch = calculate_pol2_metrics.out.pausing

  // ───────── 3.12b Aggregate list for cohort ─────────────────────────────
  def samples_lines_raw = pol2_gene_ch.map { sid, genes, c, t, r ->
    def p = file(genes)
    tuple(sid, c ?: 'NA', t ?: 'NA', (r ?: '1') as String, p.exists() ? p.toString() : null)
  }

  def samples_lines_ok = samples_lines_raw
    .filter { sid, c, t, r, path ->
      if (!path) { log.warn "[pol2-aggregate] skipping: ${sid} — missing per-sample TSV path"; return false }
      true
    }
    .map { sid, c, t, r, path -> "${sid}\t${c}\t${t}\t${r}\t${path}" }

  def samples_tsv_raw =
    Channel.of('sample_id\tcondition\ttimepoint\treplicate\tfile')
      .concat(samples_lines_ok)
      .collectFile(
        name: 'samples.tsv',
        storeDir: "${params.output_dir}/09_pol2_aggregate",
        newLine: true
      )

  def samples_tsv =
    samples_tsv_raw.map { p ->
      def f = file(p)
      def txt = f.text.replace('\uFEFF','').replace('\r','')
      def lines = txt.readLines()
      if (!lines) throw new RuntimeException("[pol2-aggregate] empty samples.tsv")

      final String HEADER = 'sample_id\tcondition\ttimepoint\treplicate\tfile'
      def out = new StringBuilder(HEADER).append('\n')

      lines.drop(1).eachWithIndex { line, idx0 ->
        if (!line?.trim()) return
        if (line.replace('\r','') == HEADER) return
        def cols = line.split('\t', -1)
        if (cols.size() != 5)
          throw new RuntimeException("[pol2-aggregate] row ${idx0+2} has ${cols.size()} fields: ${line}")
        if (!cols[4] || cols[4] == 'file')
          throw new RuntimeException("[pol2-aggregate] row ${idx0+2} has invalid file path: ${cols[4]}")
        out.append(line).append('\n')
      }

      f.text = out.toString()
      if (out.toString().trim() == HEADER)
        throw new RuntimeException("[pol2-aggregate] no valid sample rows after sanitization")
      f
    }

  summarize_pol2_metrics( samples_tsv )

  // ───────── 3.13 QC ─────────────────
  qc_pol2_tracktx(
    aligned_ch.map { sid, bam, _all, _spike, c, t, r -> tuple(sid, bam, c, t, r) }
  )
  def qc_json_meta_ch = qc_pol2_tracktx.out.json_meta

  // ───────── 3.14 REPORT BUNDLE (deterministic per-SID joins) ────────────
  def meta_ch = aligned_ch.map { sid, _fb, _ab, _sb, c, t, r -> tuple(sid, c ?: 'NA', t ?: 'NA', r ?: '1') }.distinct()

  def meta_kv      = meta_ch.map { sid, c,t,r -> tuple(sid, tuple(c,t,r)) }
  def div_kv_rep   = divergent_tx_ch          .map { sid, p, c,t,r -> tuple(sid, p) }
  def fsum_kv_rep  = functional_regions_sum_ch.map { sid, p, c,t,r -> tuple(sid, p) }
  def dens_kv_rep  = pol2_density_ch          .map { sid, p, c,t,r -> tuple(sid, p) }
  def paus_kv_rep  = pol2_pausing_ch          .map { sid, p, c,t,r -> tuple(sid, p) }
  def norm_kv_rep  = norm_factors_ch          .map { sid, p, c,t,r -> tuple(sid, p) }
  def dedu_kv_rep  = dedup_stats_ch           .map { sid, p, c,t,r -> tuple(sid, p) }
  def qcjs_kv_rep  = qc_json_meta_ch          .map { sid, p, c,t,r -> tuple(sid, p) }

  // helper: build a path under results and return string if it exists, else ''
  def maybePath = { String s ->
    def f = file(s)
    f.exists() ? f.toString() : ''
  }

  // Core report payload as KV (keep Files for path(...) inputs)
  def rep_core_kv =
      meta_kv.join(div_kv_rep)
             .join(fsum_kv_rep)
             .join(dens_kv_rep)
             .join(paus_kv_rep)
             .join(norm_kv_rep)
             .join(dedu_kv_rep)
             .join(qcjs_kv_rep)
             .map { sid, M, div, fsum, dens, paus, norm, dedu, qcj ->
               def (cond,tp,rep) = M
               tuple(sid, tuple(file(div), file(fsum), file(dens), file(paus), file(norm), file(dedu), file(qcj), cond, tp, rep))
             }

  // Normalized BWs + RAW allMap 3′ bedGraphs for the HTML extras
  def extras_kv = meta_ch.map { sid, _c,_t,_r ->
    def normDir   = "${params.output_dir}/05_normalized_tracks/${sid}"
    def bw_pos3   = maybePath("${normDir}/3p/main/cpm/pos.bw")
    def bw_neg3   = maybePath("${normDir}/3p/main/cpm/neg.bw")
    def bw_all_p3 = maybePath("${normDir}/3p/allMap/cpm/pos.bw")
    def bw_all_n3 = maybePath("${normDir}/3p/allMap/cpm/neg.bw")

    def tracksDir = "${params.output_dir}/03_genome_tracks/${sid}"
    def raw_ap3   = maybePath("${tracksDir}/3p/${sid}.allMap.3p.pos.bedgraph")
    def raw_an3   = maybePath("${tracksDir}/3p/${sid}.allMap.3p.neg.bedgraph")

    tuple(sid, tuple(raw_ap3, raw_an3, bw_pos3, bw_neg3, bw_all_p3, bw_all_n3))
  }

  // Flatten to the tuple generate_reports expects
  def rep_core_ext =
      rep_core_kv.join(extras_kv)
        .map { sid, CORE, EXS ->
          def (div, fsum, dens, paus, norm, dedu, qcj, cond, tp, rep) = CORE
          def (raw_ap3, raw_an3, BW_P3, BW_N3, BW_AP3, BW_AN3) = EXS
          tuple(sid, div, fsum, dens, paus, norm, dedu, qcj,
                     raw_ap3, raw_an3, BW_P3, BW_N3, BW_AP3, BW_AN3,
                     cond, tp, rep)
        }

  if (params.debug_report) {
    rep_core_ext.view { sid, div, fsum, dens, paus, norm, dedu, qc, ap, an, bwp, bwn, bwap, bwan, c,t,r ->
      "DEBUG[report-in] ${sid} div=${div} fsum=${fsum} dens=${dens} paus=${paus} norm=${norm} dedu=${dedu} qc=${qc} ap=${ap} an=${an} pos3bw=${bwp} neg3bw=${bwn} allMapP=${bwap} allMapN=${bwan}"
    }
  }

  // ───────── 3.15 Reports ───────────
  generate_reports( rep_core_ext )

  def report_json_ch     = generate_reports.out[2]
  def per_sample_inputs  = report_json_ch.toList()
  combine_reports( per_sample_inputs )

}
