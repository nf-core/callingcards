/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCallingcards.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome'
include { PREPARE_READS      } from '../subworkflows/local/yeast/prepare_reads'
include { ALIGN              } from '../subworkflows/local/align'
include { PROCESS_ALIGNMENTS } from '../subworkflows/local/yeast/process_alignments'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if(params.fasta){
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).collect()
} else {
    exit 1, 'Either a valid configured `genome` or a `fasta` file must be specified.'
}

if(params.gtf){
    ch_gtf = Channel.fromPath(params.gtf, checkIfExists:true).collect()
} else {
    exit 1, 'Either a valid configured `genome` or a `gtf` file must be specified.'
}

ch_regions_mask = params.regions_mask ?
        Channel.fromPath(params.regions_mask, checkIfExists: true)
                .collect().map{ it -> [[id:it[0].getSimpleName()], it[0]]} :
        Channel.empty()

additional_fasta = params.additional_fasta ?
        Channel.fromPath(params.additional_fasta, checkIfExists: true).collect() :
        Channel.empty()

def rseqc_modules = params.rseqc_modules ?
    params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } :
    []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CALLINGCARDS_YEAST {

    // instantiate channels
    ch_versions          = Channel.empty()
    ch_fasta_index       = Channel.empty()
    ch_bam_index         = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flatstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()

    //
    // VALIDATE AND PARSE INPUTS
    //
    Channel.fromSamplesheet("input")
        .multiMap{ sample, fastq_1, fastq_2, barcode_details ->
            def single_end = fastq_2.size() == 0
            if (single_end){
                exit "Currently the yeast pipeline expects " +
                "paired end reads due to the barcodes."
            }
            def meta = ["id": sample, "single_end": single_end]
            reads: [meta, [fastq_1, fastq_2]]
            barcode_details: [meta, barcode_details]}
        .set{ ch_input }

    // SUBWORKFLOW_2: Index the genome in various ways
    //
    // input: fasta file (genome), a regions mask (bed format), additional
    //        sequences (fasta format) to append to the genome after masking,
    //        and a gtf file
    // output: fasta (masked fasta), fai (fasta index),
    //         genome_bed (gtf in bed format), bwamem2_index (bwa mem2 index),
    //         bwa_index (bwa aln index), bowtie_index (bowtie index),
    //         bowtie2_index (bowtie2 index), versions
    //         NOTE: the aligner index channels will be empty except for the
    //         aligner specified in params.aligner
    PREPARE_GENOME(
        ch_fasta,
        ch_regions_mask,
        additional_fasta,
        ch_gtf
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW_3: run sequencer level QC, extract barcodes and trim
    // input: [ val(meta), [ path(fastq1), path(fastq2) ] ] where the reads
    //        list may be [ path(fastq1) ] if the input is single-end
    // output: 'reads' with structure [ val(meta),
    //                                  [ path(fastq1), path(fastq2) ] ]
    //         '(raw/demux)_fastqc_html' with structure [ val(meta),
    //                                                    path(fastqc_html) ]
    //         '(raw/demux)_fastqc_zip' with structure [ val(meta),
    //                                                   path(fastqc_zip) ]
    //         'trimmomatic_log' with structure [ val(meta),
    //                                            path(trimmomatic_log) ]
    //         'versions' with structure [ path(versions) ]
    //
    PREPARE_READS (
        ch_input.reads,
        ch_input.barcode_details
    )
    ch_versions = ch_versions.mix(PREPARE_READS.out.versions)

    //
    // SUBWORKFLOW_4: align reads
    // input: post-processed reads. note that the fastq have been split into
    //        chunks based on params.split_fastq_chunk_[size/part]. the metadata
    //        includes a key value pair split: <split_number>, eg split: 1
    // output: channel 'bam' with structure [ val(meta), path(bam), path(bai) ]
    //         channel 'versions' with structure [ path(versions) ]
    ALIGN (
        PREPARE_READS.out.reads,
        PREPARE_GENOME.out.bwamem2_index,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.bowtie_index
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    // create a channel from the ch_input barcode details
    // of structure [ [id: sample_id], path(barcode_details) ]
    // this will be used to join the barcode_details to the alignments
    ch_input.barcode_details
        .map{ meta, barcode_details ->
                [meta.id, barcode_details]}
        .set{ ch_barcode_details_join }

    // join the alignments to the barcode details on the shared sample id
    ALIGN.out.bam
        .map{ meta, bam, bai ->
                [meta.id, meta, bam, bai] }
        .combine(ch_barcode_details_join, by:0)
        .map{ key, meta, bam, bai, barcode_details ->
                [meta, bam, bai, barcode_details] }
        .set{ ch_aln_with_details }

    //
    // SUBWORKFLOW_5: process the alignnment into a qbed file
    // input: 'bam' with structure [ val(meta), path(bam),
    //                               path(bai), path(barcode_details) ],
    //       'fasta' (genome fasta), 'genome_bed' (genome bed),
    //       'rseqc_modules' (list of rseqc modules to run),
    //       'fai' [val(meta), path(fai)], 'gtf' (gtf file),
    //       'biotypes_header_multiqc' (biotypes header for multiqc)
    // output: All QC outputs have structure [ val(meta), path(file) ]
    //         'samtools_stats', 'samtools_flatstat', 'samtools_idxstats',
    //         'picard_qc', 'rseqc_bamstat', 'rseqc_infer_experiment',
    //         'rseqc_inner_distance', 'rseqc_read_distribution',
    //         'rseqc_read_duplication', 'rseqc_tin',
    //         'versions'
    PROCESS_ALIGNMENTS (
        ch_aln_with_details,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.genome_bed,
        rseqc_modules,
        ch_gtf
    )
    ch_versions = ch_versions.mix(PROCESS_ALIGNMENTS.out.versions)


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCallingcards.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCallingcards.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.raw_fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.demux_fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.trimmomatic_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.samtools_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.samtools_flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.samtools_idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.picard_qc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.rseqc_bamstat.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.rseqc_inferexperiment.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.rseqc_innerdistance.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.rseqc_readdistribution.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.rseqc_readduplication.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(PROCESS_ALIGNMENTS.out.rseqc_tin.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
