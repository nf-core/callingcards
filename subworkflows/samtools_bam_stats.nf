//
// Run SAMtools stats, flagstat and idxstats
// COPIED FROM NF-CO/RNASEQ
//

include { SAMTOOLS_STATS    } from "${projectDir}/modules/nf-core/samtools/stats/main"
include { SAMTOOLS_IDXSTATS } from "${projectDir}/modules/nf-core/samtools/idxstats/main"
include { SAMTOOLS_FLAGSTAT } from "${projectDir}/modules/nf-core/samtools/flagstat/main"

workflow SAMTOOLS_BAM_STATS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai/csi] ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_STATS ( ch_bam_bai, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
