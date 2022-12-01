//
// Read QC, UMI extraction and trimming
// NOTE: copied from nf-co/rnaseq/subworkflows/nf-core/fastqc_umitools_trimgalore.nf
//

// TODO this is the longest process -- split into n groups and do in parallel
// use this as template https://www.nextflow.io/example3.html

include { FASTQC           } from "${projectDir}/modules/nf-core/fastqc/main"
include { UMITOOLS_EXTRACT } from "${projectDir}/modules/nf-core/umitools/extract/main"
include { TRIMMOMATIC      } from "${projectDir}/modules/nf-core/trimmomatic/main"

workflow PREPARE_READS {
    take:
    reads         // channel: [ val(meta), [ reads ] ]

    main:

    // log output
    ch_versions = Channel.empty()
    umi_log   = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    trimmomatic_log = Channel.empty()

    // read output
    ch_reads = Channel.empty()


    // run umi extract to add barcodes to fastq id lines
    UMITOOLS_EXTRACT ( reads )
    umi_log     = UMITOOLS_EXTRACT.out.log
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)

    // run fastqc after trimming off the barcodes, etc
    FASTQC ( reads ).html.set { fastqc_html }
    fastqc_zip  = FASTQC.out.zip
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // if the reduce_to_se flag is set, the user wishes to reduce paired
    // end reads to R1 only
    if(params.reduce_to_se){
        // reduce the reads from umitools extract to
        // [ val(meta), path(...R1.fastq.gz)]
        // where the meta.single_end tag is set to true
        UMITOOLS_EXTRACT.out.reads
            .map{ meta, reads ->
                [to_single_end(meta), reads[0]]}
            .set{ trimmomatic_input }
        // if the user also wishes to crop r1 (eg, reduce to 88 bp at most)
        // then do so
        if(params.r1_crop){
            TRIMMOMATIC(
                trimmomatic_input
            )
            ch_reads = ch_reads.mix(TRIMMOMATIC.out.trimmed_reads)
            trimmomatic_log = trimmomatic_log.mix(TRIMMOMATIC.out.log)
            ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
        // else, return the uncropped SE read data to ch_reads
        }else {
            ch_reads = ch_reads.mix(trimmomatic_input)
        }
    // if reduce_to_se is NOT set, still allow r1_crop
    // TODO this is extremely error prone and needs to be fixed -- need to be
    // able to check if meta.single_end is set
    } else if (params.r1_crop){
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
    reads = ch_reads       // channel: [ val(meta), [ reads ] ]
    fastqc_html            // channel: [ val(meta), [ html ] ]
    fastqc_zip             // channel: [ val(meta), [ zip ] ]
    umi_log                // channel: [ val(meta), [ log ] ]
    trimmomatic_log        // channel: [ val(meta), [ log ] ]
    versions = ch_versions // channel: [ versions.yml ]
}

def to_single_end(Map meta) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.single_end = true

    return new_meta
}
