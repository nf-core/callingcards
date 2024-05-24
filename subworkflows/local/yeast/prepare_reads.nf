//
// Read QC, UMI extraction and trimming
// NOTE: inspired by nf-co/rnaseq/subworkflows/nf-core/fastqc_umitools_trimgalore.nf
//

include { FASTQC as FASTQCRAW   } from "../../../modules/nf-core/fastqc/main"
include { FASTQC as FASTQCDEMUX } from "../../../modules/nf-core/fastqc/main"
include { SEQKIT_SPLIT2         } from "../../../modules/nf-core/seqkit/split2/main"
include { TRIMMOMATIC           } from "../../../modules/nf-core/trimmomatic/main"
include { DEMULTIPLEX           } from "../../../modules/local/callingcardstools/yeast/demultiplex/main"
include { CONCATQC              } from "../../../modules/local/callingcardstools/yeast/concatQC/main"
include { CONCATFASTQ           } from "../../../modules/local/concatFastq/main"

workflow PREPARE_READS {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    barcode_details // channel: [ val(meta), [ barcode_details ] ]

    main:

    // log output
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    trimmomatic_log = Channel.empty()

    // read output
    demux_reads = Channel.empty()
    output_reads = Channel.empty()

    // run FASTQC on the raw input reads
    FASTQCRAW ( reads )
    // NOTE: version is stored in the demultiplexed reads fastqc step below

    FASTQCRAW.out.html.set{ raw_fastqc_html }
    FASTQCRAW.out.zip.set{ raw_fastqc_zip }

    SEQKIT_SPLIT2 ( reads )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    // NOTE!! I am relying on groupTuple() to put the reads in the right
    // order. need to do thi explicitly somehow.
    SEQKIT_SPLIT2.out.reads
        .transpose()
        .map{ meta, read ->
            [add_split(meta, read.getName()), read]}
        .groupTuple()
        .map{meta, reads ->
                [meta.id, meta, reads]}
        .combine(barcode_details.map{meta, barcode_json -> [meta.id, barcode_json]}, by:0)
        .map{id, meta, read, barcode_details ->
            [meta, read, barcode_details]}
        .set{ split_reads_with_barcode_details }

    DEMULTIPLEX(
        split_reads_with_barcode_details
    )
    ch_versions = ch_versions.mix(DEMULTIPLEX.out.versions)

    // transpose associates the metadata with each individual read
    // the mapping adds the tf info, removes the split info, and adds
    // the read end
    // groupTuple groups the reads by metadata, so all splits of R1 for a given
    // tf in a given batch will be grouped together
    // next we map to and sort the reads by the split number to ensure the
    // order is correct
    DEMULTIPLEX.out.reads
        .transpose()
        .map{meta, read ->
            [add_tf_remove_split_add_read(meta, read.getName()), read]}
        .groupTuple()
        .map{ meta, reads ->
                [meta, sortFilesBySplit(reads)]}
        .set{ demux_reads }

    // concatenate the splits by TF and read end
    CONCATFASTQ(
        demux_reads
    )
    ch_versions = ch_versions.mix(CONCATFASTQ.out.versions)

    // group all split by id and TF to collect the barcode QC picke objects to
    // ch_demux_qc, which will be input into the CONCATQC step
    DEMULTIPLEX.out.pickle
        .map{meta, pickle ->
            [[id:meta.id, single_end: meta.single_end], pickle]}
        .groupTuple()
        .join(barcode_details, by:0)
        .map{meta, pickle, barcode_details ->
            [meta, pickle, barcode_details]}
        .set{ ch_demux_qc }

    // combine the split barcode QC objects into a single object for each
    // sample_id + TF
    CONCATQC(
        ch_demux_qc
    )
    ch_versions = ch_versions.mix(CONCATQC.out.versions)

    // first remove the read end info from the metadata,
    // next groupTuple to create a [metadata, [read1, read2]] entry by
    // TF. Note that the reads are sorted so that R1 preceeds R2
    CONCATFASTQ.out.reads
        .map{
            meta, reads ->
                [remove_read_end(meta), reads]
        }
        .groupTuple()
        .map{meta, reads ->
                [meta, sortFilesByReadEnd(reads)]}
        .filter{meta, reads ->
            reads.every {it.size() > 0} }
        .set{ demux_concat_reads }

    // switch single_end to true in meta and select only R1
    demux_concat_reads
        .map{ meta, reads ->
            [to_single_end(meta), reads[0]]}
            .set{ trimmomatic_input }

    // trim the end of the reads based on params.r1_crop
    TRIMMOMATIC( trimmomatic_input  )
    trimmomatic_log = trimmomatic_log.mix(TRIMMOMATIC.out.log)
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

    output_reads = output_reads.mix(TRIMMOMATIC.out.trimmed_reads)

    // run fastqc by TF after trimming (note that the barcodes are
    // removed at this point )
    FASTQCDEMUX ( output_reads )
    FASTQCDEMUX.out.html.set{ demux_fastqc_html }
    FASTQCDEMUX.out.zip.set { demux_fastqc_zip  }
    ch_versions = ch_versions.mix(FASTQCDEMUX.out.versions.first())

    emit:
    reads = output_reads   // channel: [ val(meta), [ reads ] ]
    raw_fastqc_html        // channel: [ val(meta), [ html ] ]
    raw_fastqc_zip         // channel: [ val(meta), [ zip ] ]
    demux_fastqc_html      // channel: [ val(meta), [ html ] ]
    demux_fastqc_zip       // channel: [ val(meta), [ zip ] ]
    trimmomatic_log        // channel: [ val(meta), [ log ] ]
    versions = ch_versions // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Groovy Utility Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
import nextflow.util.ArrayBag

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

def to_single_end(Map meta) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.single_end = true

    return new_meta
}

// Extract R1 or R2 from the file name
def extractReadType(String path) {
    // Regex pattern to match R1 or R2 in the file name
    def pattern = /(R1|R2)/

    // Extract R1 or R2
    def matcher = path =~ pattern
    if (matcher.find()) {
        return matcher[0][0]
    } else {
        return null
    }
}

// add the tf info to the meta
def add_tf_remove_split_add_read(Map meta, String read) {
    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.remove('split')

    String tf = (read =~ /^(.+?)_/)[0][1]

    new_meta.tf = tf

    new_meta.read_end = extractReadType(read)

    return new_meta
}

// remove read_end
def remove_read_end(Map meta) {

    meta.remove('read_end')

    return meta
}

// sort split files by the split digit in the iflename
def sortFilesBySplit(ArrayBag<Path> files) {
    // Define a closure that extracts the digit surrounded by underscores
    def extractDigit = { String path ->
        def pattern = /_(\d+)_/
        def matcher = path =~ pattern
        if (matcher.find()) {
            return matcher[0][1].toInteger()
        } else {
            return null
        }
    }

    // Sort the ArrayBag using the closure
    files.sort { a, b ->
        extractDigit(a.toString()) <=> extractDigit(b.toString())
    }
}

def sortFilesByReadEnd(ArrayBag<Path> files) {
    // Define a closure that extracts the read end (R1 or R2) from the file path
    def extractReadEnd = { String path ->
        def pattern = /_(R[12])_/
        def matcher = path =~ pattern
        if (matcher.find()) {
            return matcher[0][1]
        } else {
            return null
        }
    }

    // Sort the ArrayBag using the closure
    files.sort { a, b ->
        extractReadEnd(a.toString()) <=> extractReadEnd(b.toString())
    }
}
