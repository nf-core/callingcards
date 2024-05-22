//
// Check input samplesheet and get read channels
//
include { BAM_STATS_SAMTOOLS            } from "../../nf-core/bam_stats_samtools/main"
include { BAM_RSEQC                     } from "../../nf-core/bam_rseqc/main"
include { COUNT_HOPS                    } from "../../../modules/local/callingcardstools/mammals/count_hops/main"
include { CONCATQC                      } from "../../../modules/local/callingcardstools/mammals/concatQC/main"
include { PICARD_COLLECTMULTIPLEMETRICS } from "../../../modules/nf-core/picard/collectmultiplemetrics/main"
include { SAMTOOLS_SORT                 } from "../../../modules/nf-core/samtools/sort/main"
include { SAMTOOLS_INDEX                } from "../../../modules/nf-core/samtools/index/main"
include { SAMTOOLS_MERGE                } from "../../../modules/nf-core/samtools/merge/main"

workflow PROCESS_ALIGNMENTS {
    take:
    // [ val(meta), path(bam), path(bai), path(barcode_details) ]
    // note that the fastq were split into chunks of equal size. meta has a
    // key: value pair split: <split_number>, eg split: 1
    aln
    fasta                   // [val(meta), path(genome.fasta)]
    fai                     // [ val(meta), path(genome.fasta.fai) ]
    ch_genome_bed           // path(genome.bed)
    rseqc_modules           // [ val(module_list) ]
    ch_gtf                  // [ path(gtf) ]

    main:

    ch_versions = Channel.empty()

    // parition alignments into passing/failing bam files and pass on the
    // pickled qc qbed and barcode qc
    COUNT_HOPS (
        aln,
        fasta.map{meta, fasta -> fasta},
        fai.map{meta,fai -> fai} )

    ch_versions = ch_versions.mix(COUNT_HOPS.out.versions)

    // combine the qbed pickles from the splits
    COUNT_HOPS.out.qbed_pkl
        .map{meta, pkl ->
        [[id:meta.id], pkl]}
        .groupTuple()
        .set{ ch_grouped_pkl }

    // combine the barcode qc pickles from the splits
    COUNT_HOPS.out.barcode_qc_pkl
        .map{meta, pkl ->
        [[id:meta.id], pkl]}
        .groupTuple()
        .set{ ch_grouped_barcode_qc_pkl }

    // combine the qbed and barcode qc pickles from the splits
    // publish the combined QC stats from this step
    CONCATQC ( ch_grouped_pkl, ch_grouped_barcode_qc_pkl )
    ch_versions = ch_versions.mix(CONCATQC.out.versions)

    // combine the passing and failing bam splits to ch_bam. To
    // each of the meta maps, add the key 'qc_status'
    // with value 'passing'/'failing'
    COUNT_HOPS.out.passing_bam
        .map{meta, bam ->
        [[id:meta.id, qc_status: 'passing'], bam]}
        .groupTuple()
        .mix(COUNT_HOPS.out.failing_bam
        .map{meta, bam ->
        [[id:meta.id, qc_status: 'failing'], bam]}
        .groupTuple())
        .set{ ch_bam }

    // merge the passing and failing sets of bams into a single bam
    SAMTOOLS_MERGE( ch_bam, fasta, fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // sort the merged bams
    SAMTOOLS_SORT( SAMTOOLS_MERGE.out.bam, [[], []] )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // index the merged, sorted bams
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // join the merged, sorted bam channel with the index channel to create
    // ch_bam_bai with structure [ val(meta), path(bam), path(bai) ]
    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ ch_bam_bai }

    // run samtools stats, flagstat and idxstats on the merged, sorted bams
    BAM_STATS_SAMTOOLS(
        ch_bam_bai,
        fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    // run picard collectmultiplemetrics on the merged, sorted bams
    PICARD_COLLECTMULTIPLEMETRICS(
        ch_bam_bai,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    // run rseqc on the merged, sorted bams
    BAM_RSEQC(
        ch_bam_bai,
        ch_genome_bed,
        rseqc_modules
    )

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

