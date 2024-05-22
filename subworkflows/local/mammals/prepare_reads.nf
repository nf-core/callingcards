//
// Read QC, UMI extraction and trimming
// NOTE: copied from nf-co/rnaseq/subworkflows/nf-core/fastqc_umitools_trimgalore.nf
//

include { FASTQC           } from "../../../modules/nf-core/fastqc/main"
include { SEQKIT_SPLIT2    } from "../../../modules/nf-core/seqkit/split2/main"
include { UMITOOLS_EXTRACT } from "../../../modules/nf-core/umitools/extract/main"
include { TRIMMOMATIC      } from "../../../modules/nf-core/trimmomatic/main"

workflow PREPARE_READS {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:

    // log output
    ch_versions     = Channel.empty()
    umi_log         = Channel.empty()
    fastqc_html     = Channel.empty()
    fastqc_zip      = Channel.empty()
    trimmomatic_log = Channel.empty()

    // run fastqc after trimming off the barcodes, etc
    FASTQC ( reads ).html.set { fastqc_html }
    fastqc_zip  = FASTQC.out.zip
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // read output
    ch_reads = Channel.empty()

    SEQKIT_SPLIT2 ( reads )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    SEQKIT_SPLIT2.out.reads
        .transpose()
        .map{ meta, read1 ->
            if (!meta.single_end){
                exit 1, "Paired-end reads are not supported for mammals data. " +
                "If you must submit both ends for barcode purposes, you will " +
                "need to open a feature request on GitHub or post on the " +
                "nf-core/callingcards slack channel. Otherwise, re-submit " +
                "with only fastq_1 in the samplesheet."
            }
            [add_split(meta, read1.getName()), [read1]]}
        .set{ ch_split_reads }

    // run umi extract to add barcodes to fastq id lines
    UMITOOLS_EXTRACT ( ch_split_reads )
    umi_log     = UMITOOLS_EXTRACT.out.log
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)

    if (params.r1_crop){
        TRIMMOMATIC(
                UMITOOLS_EXTRACT.out.reads
            )
            ch_reads = ch_reads.mix(TRIMMOMATIC.out.trimmed_reads)
            trimmomatic_log = trimmomatic_log.mix(TRIMMOMATIC.out.log)
            ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    } else {
        ch_reads = ch_reads.mix(UMITOOLS_EXTRACT.out.reads)
    }

    emit:
    // this is performed prior to any umi extraction or trimming
    fastqc_html            // channel: [ val(meta), [ html ] ]
    fastqc_zip             // channel: [ val(meta), [ zip ] ]
    // note that the meta now includes the key split and value split_number (eg split: 1)
    reads = ch_reads       // channel: [ val(meta), [ reads ] ]
    umi_log                // channel: [ val(meta), [ log ] ]
    trimmomatic_log        // channel: [ val(meta), [ log ] ]
    // software versions
    versions = ch_versions // channel: [ versions.yml ]
}

// Groovy functions

def extractDigitBeforeExtension(String path) {
    // Regex pattern to match the digit before the file extension
    def pattern = /(\d+)(?=\.(fastq|fastq\.gz|fq|fq\.gz)$)/

    // Extract the digit
    def matcher = path =~ pattern
    if (matcher.find()) {
        return matcher[0][1].toInteger()
    } else {
        return null
    }
}

def add_split(Map meta, String read){
    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.split = extractDigitBeforeExtension(read)

    return new_meta
}
