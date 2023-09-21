//
// Check input samplesheet and get read channels
//
include { BAM_STATS_SAMTOOLS            } from "${projectDir}/subworkflows/nf-core/bam_stats_samtools/main"
include { PICARD_COLLECTMULTIPLEMETRICS } from "${projectDir}/modules/nf-core/picard/collectmultiplemetrics/main"
include { BAM_RSEQC                     } from "${projectDir}/subworkflows/nf-core/bam_rseqc/main"
include { HOPS_MAMMALS                  } from "${projectDir}/modules/local/hops_mammals/main"
include { SAMTOOLS_INDEX                } from "${projectDir}/modules/nf-core/samtools/index/main"

workflow PROCESS_ALIGNMENTS {
    take:
    aln //channel: [ val(meta), path(bam), path(bai), path(barcode_details) ]
    fasta // path(genome.fasta)
    ch_genome_bed
    rseqc_modules
    fai // [ val(meta), path(genome.fasta.fai) ]
    ch_gtf // [ path(gtf) ]

    main:

    ch_versions = Channel.empty()

    // create channel without bai for those steps that take it as input

    aln
      .map{meta,bam,bai,barcode_details -> [meta,bam,bai] }
      .set{ bam_bai }

    BAM_STATS_SAMTOOLS(
        bam_bai,
        fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    PICARD_COLLECTMULTIPLEMETRICS(
        bam_bai,
        fasta.map{it -> ['',it]},
        fai
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    BAM_RSEQC(
        bam_bai,
        ch_genome_bed,
        rseqc_modules
    )

    HOPS_MAMMALS (
        aln,
        fasta,
        fai.map{meta,fai -> fai}
    )
    ch_versions = ch_versions.mix(HOPS_MAMMALS.out.versions)

    HOPS_MAMMALS.out.passing_bam
        .mix(HOPS_MAMMALS.out.failing_bam)
        .set{ ch_parse_bam_out }

    SAMTOOLS_INDEX(
        ch_parse_bam_out
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    aln
    .map{ meta,bam,bai,barcode -> [meta,bam] }
    .combine(ch_gtf)
    .set{ ch_featurecounts_input }

    emit:
    samtools_stats = BAM_STATS_SAMTOOLS.out.stats
    samtools_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    samtools_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
    picard_qc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    rseqc_bamstat            = BAM_RSEQC.out.bamstat_txt
    rseqc_inferexperiment    = BAM_RSEQC.out.inferexperiment_txt
    rseqc_innerdistance      = BAM_RSEQC.out.innerdistance_freq
    rseqc_junctionannotation = BAM_RSEQC.out.junctionannotation_log
    rseqc_junctionsaturation = BAM_RSEQC.out.junctionsaturation_rscript
    rseqc_readdistribution   = BAM_RSEQC.out.readdistribution_txt
    rseqc_readduplication    = BAM_RSEQC.out.readduplication_pos_xls
    rseqc_tin                = BAM_RSEQC.out.tin_txt
    versions = ch_versions
}
