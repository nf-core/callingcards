//
// Process alignments
//
include { BAM_STATS_SAMTOOLS            } from "../../nf-core/bam_stats_samtools/main"
include { BAM_RSEQC                     } from "../../nf-core/bam_rseqc/main"
include { COUNT_HOPS                    } from "../../../modules/local/callingcardstools/yeast/count_hops/main"
include { PICARD_COLLECTMULTIPLEMETRICS } from "../../../modules/nf-core/picard/collectmultiplemetrics/main"
include { SAMTOOLS_INDEX                } from "../../../modules/nf-core/samtools/index/main"

workflow PROCESS_ALIGNMENTS {
    take:
    aln            // [ val(meta), path(bam), path(bai), path(barcode_details) ]
    fasta          // [val(meta), path(genome.fasta)]
    fai            // [ val(meta), path(genome.fasta.fai) ]
    ch_genome_bed  // [ path(genome.bed) ]
    rseqc_modules  // [ val(module_list) ]
    ch_gtf         // [ path(gtf) ]

    main:

    ch_versions = Channel.empty()

    // split the alignment file into passing/failing. output alignment level
    // qc metrics
    COUNT_HOPS (
        aln,
        fasta.map{meta, fasta -> fasta},
        fai.map{meta,fai -> fai}
    )
    ch_versions = ch_versions.mix(COUNT_HOPS.out.versions)

    // add the qc_status key to the passing and failing hops channels
    // respectively. Note that input should be sorted, and output will
    // retain sorted order
    COUNT_HOPS.out.passing_bam
        .map{meta, bam, bai ->
            [meta, add_qc_status(meta, 'passing'), bam, bai]}
        .set{ ch_passing_bam }

    COUNT_HOPS.out.failing_bam
        .map{meta, bam, bai ->
            [meta, add_qc_status(meta, 'failing'), bam, bai]}
        .set{ ch_failing_bam }


    // mix the passing and failing channels together to perform QC
    // NOTE: the filter clause removes from the channel any items with
    // no reads in the bam file
    ch_passing_bam
        .mix(ch_failing_bam)
        .map{meta_join, meta, bam, bai ->
            [meta, bam, bai]}
        .set{ ch_parse_bam_out }

    // run samtools stats, flagstat and idxstats on the merged, sorted bams
    BAM_STATS_SAMTOOLS(
        ch_parse_bam_out,
        fasta,
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    // run picard collectmultiplemetrics on the merged, sorted bams
    PICARD_COLLECTMULTIPLEMETRICS(
        ch_parse_bam_out,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    // run rseqc on the merged, sorted bams
    BAM_RSEQC(
        ch_parse_bam_out,
        ch_genome_bed,
        rseqc_modules
    )
    ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

    // TODO: add featurecounts
    // aln
    // .map{ meta,bam,bai,barcode -> [meta,bam] }
    // .combine(ch_gtf)
    // .set{ ch_featurecounts_input }

    emit:
    samtools_stats           = BAM_STATS_SAMTOOLS.out.stats
    samtools_flagstat        = BAM_STATS_SAMTOOLS.out.flagstat
    samtools_idxstats        = BAM_STATS_SAMTOOLS.out.idxstats
    picard_qc                = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    rseqc_bamstat            = BAM_RSEQC.out.bamstat_txt
    rseqc_inferexperiment    = BAM_RSEQC.out.inferexperiment_txt
    rseqc_innerdistance      = BAM_RSEQC.out.innerdistance_freq
    rseqc_readdistribution   = BAM_RSEQC.out.readdistribution_txt
    rseqc_readduplication    = BAM_RSEQC.out.readduplication_pos_xls
    rseqc_tin                = BAM_RSEQC.out.tin_txt
    versions                 = ch_versions
}

// add the split number to the metadata; keep all other key:value pairs in meta
def add_qc_status(Map meta, String status){
    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.qc_status = status

    return new_meta
}
