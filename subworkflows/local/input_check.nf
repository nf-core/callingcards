//
// Check input samplesheet and get read channels
//
include { SAMPLESHEET_CHECK } from "../../modules/local/samplesheet_check"

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .multiMap { it ->
            reads: create_fastq_channel(it, params.datatype)
            barcode_details: create_barcode_details_channel(it, params.datatype) }
        .set { result }

    emit:
    // NOTE: there is a subtlety in how single_end is set. If the
    // datatype is mammals, then single_end is set to true for all samples
    // and a warning is written to std out for the user to the effect that
    // the mammals pipeline only processes R1 reads at this point.
    reads           = result.reads                    // channel: [ val(meta), [ reads ] ]
    barcode_details = result.barcode_details          // channel: [ val(meta), [barcode_details]]
    versions        = SAMPLESHEET_CHECK.out.versions  // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row, String datatype) {
    // create meta map
    def meta        = [:]
    meta.id         = row.sample
    meta.single_end = !row.single_end.toBoolean() & datatype=='mammal' ?
        true :
        row.single_end.toBoolean()

    if(datatype == 'mammal' & !row.single_end.toBoolean()){
        log.warn"Removing R2 from mammals input ${meta.id} and setting " +
            "meta.single_end to true. The mammals workflow " +
            "does not use R2 reads and they are not currently used " +
            "for alignment."
    }

    // this is done to make the test profiles work. test data are in the
    // projectDir currently
    // def fastq_1 = file(row.fastq_1).exists() ?
    //                 row.fastq_1 :
    //                 "${projectDir}/${row.fastq_1}"
    fastq_1 = row.fastq_1

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(fastq_1) ] ]
    } else {
    // this is done to make the test profiles work. test data are in the
    // projectDir currently
        // def fastq_2 = file(row.fastq_2).exists() ?
        //         row.fastq_2 :
        //         "${projectDir}/${row.fastq_2}"
        fastq_2 = row.fastq_2

        if (!file(fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${fastq_2}"
        }

        fastq_meta = [ meta, [ file(fastq_1), file(fastq_2) ] ]
    }
    return fastq_meta
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_barcode_details_channel(LinkedHashMap row, String datatype) {
    // create meta map
    def meta        = [:]
    meta.id         = row.sample
    meta.single_end = !row.single_end.toBoolean() & datatype=='mammal' ?
        true :
        row.single_end.toBoolean()

    // if the row has a tf column, add it to the meta map. this will
    // mean that the demultiplexing step is skipped
    if (row.tf) {
        meta.tf = row.tf
    }

    // this is done to make the test profiles work. test data are in the
    // projectDir currently
    // def barcode_details = file(row.barcode_details).exists() ?
    //                 row.barcode_details :
    //                 "${projectDir}/${row.barcode_details}"
    barcode_details = row.barcode_details

    // add path(s) of the fastq file(s) to the meta map
    def barcode_details_meta = []
    if (!file(barcode_details).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> barcode_details file does not exist!\n${barcode_details}"
    }
    barcode_details_meta = [ meta, [ file(barcode_details) ] ]
    return barcode_details_meta
}
