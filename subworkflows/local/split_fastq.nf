//
// Check input samplesheet and get read channels
//

include { NF-CORE_TRIMMOMATIC } from '../modules/nf-core/modules/nf-core/trimmomatic/main'
include { PARSE_FASTQ } from '../../modules/local/parse_fastq'

workflow INPUT_CHECK {
    take:
    reads // [ val(meta), [ reads ] ]

    main:
    PARSE_FASTQ ( reads )

    ch_split_r1     = Channel.empty()
    trimmomatic_log = Channel.empty()
    ch_versions     = Channel.empty()

    PARSE_FASTQ.out.meta
        .combine(PARSE_FASTQ.out.reads.flatten())
        .map{ meta, r1 -> reformat_meta(meta,r1),r1}
        .set{parsed_r1}

    ch_versions = ch_versions.mix(PARSE_FASTQ.out.versions)

    if (params.r1_crop){
        TRIMMOMATIC (
            parsed_r1
        )
        ch_split_r1 = ch_split_r1.mix(TRIMMOMATIC.out.trimmed_reads)
        trimmomatic_log = trimmomatic_log.mix(TRIMMOMATIC.out.log)
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    } else {
        ch_split_r1 = ch_split_r1.mix(parsed_r1)
    }

    ch_split_r1 = ch_split_r1
        .combine(parsed_barcode_details)

    emit:
    r1 = ch_split_r1       // channel: [ val(meta), path(r1), path(barcode_details) ]
    versions = ch_versions // channel: [ versions.yml ]
}

def reformat_meta(Map meta,r1) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    def (_,split) = (r1 =~ /.*(?=_R1.fq)/ )[0]

    new_meta['single_end'] = true
    new_meta["split"] = split

    return new_meta
}
