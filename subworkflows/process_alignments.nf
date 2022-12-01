//
// Check input samplesheet and get read channels
//
include { SAMTOOLS_BAM_STATS            } from "${projectDir}/subworkflows/samtools_bam_stats"
include { PICARD_COLLECTMULTIPLEMETRICS } from "${projectDir}/modules/nf-core/picard/collectmultiplemetrics/main"
include { PRESEQ_CCURVE                 } from "${projectDir}/modules/nf-core/preseq/ccurve/main"
include { PARSEBAM                      } from "${projectDir}/modules/local/parse_bam/main"

workflow PROCESS_ALIGNMENTS {
    take:
    aln //channel: [ val(meta), path(bam), path(bai), path(barcode_details) ]
    fasta // path(genome.fasta)
    fai // [ val(meta), path(genome.fasta.fai) ]

    main:

    ch_versions = Channel.empty()

    // create channel without bai for those steps that take it as input

    aln
      .map{meta,bam,bai,barcode_details -> [meta,bam,bai] }
      .set{ bam_bai }

    SAMTOOLS_BAM_STATS(
        bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_STATS.out.versions)

    PICARD_COLLECTMULTIPLEMETRICS(
        bam_bai.map{meta,bam,bai -> [meta,bam]},
        fasta,
        fai.map{meta,fai -> fai}
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    PRESEQ_CCURVE(
        bam_bai.map{meta,bam,bai -> [meta,bam]}
    )
    ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions)

    PARSEBAM (
        aln,
        fasta,
        fai.map{meta,fai -> fai}
    )
    ch_versions = ch_versions.mix(PARSEBAM.out.versions)

    emit:
    samtools_stats = SAMTOOLS_BAM_STATS.out.stats
    samtools_flagstat = SAMTOOLS_BAM_STATS.out.flagstat
    samtools_idxstats = SAMTOOLS_BAM_STATS.out.idxstats
    picard_qc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    preseq_ccurve = PRESEQ_CCURVE.out.c_curve
    versions = ch_versions
}
