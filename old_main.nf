#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/callingcards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/callingcards
    Website: https://nf-co.re/callingcards
    Slack  : https://nfcore.slack.com/channels/callingcards
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf = WorkflowMain.getGenomeAttribute(params, 'gtf')

if (params.datatype == 'yeast'){
    if(params.genome == 'R64-1-1'){
        log.info"${projectDir}"
        params.regions_mask = "${projectDir}/assets/yeast/igenomes/R64-1-1/regions_mask.bed"
        log.info"Using default regions mask for yeast analysis: ${params.regions_mask}"

        params.fasta_index = null
        log.info"Genome fasta will be indexed in the workflow"
    }
    if(!params.containsKey('additional_fasta')){
        params.additional_fasta = "${projectDir}/assets/yeast/plasmid_sequences.fasta"
        log.info"Using default plasmid sequences for yeast analysis: ${params.additional_fasta}"
    }
}

if (params.genome == 'GRCm38'){
    if(params.aligner == 'bwa'){
        params.bwa_index = WorkflowMain.getGenomeAttribute(params, 'bwa')
    }
}

if(!params.containsKey('regions_mask')){
    params.regions_mask = null
    log.info"Regions mask not specified. The entire genome will be used for alignment"
}

if(!params.containsKey('additional_fasta')){
    params.additional_fasta = null
    log.info"Additional fasta not specified. No additional sequences beyond those in the genome fasta will be used for alignment"
}

if (!params.containsKey('fasta_index')){
    params.fasta_index = null
    log.info"Genome fasta will be indexed in the workflow"
}

// these parameters are used in the pipeline and may or may not be set either
// by the user or through the --genome argument. If they are null, then
// the appropriate index will be created in the workflow
if(!params.containsKey('bwa_index')){
    params.bwa_index = null
}
if(!params.containsKey('bwamem2_index')){
    params.bwamem2_index = null
}
if(!params.containsKey('bowtie_index')){
    params.bowtie_index = null
}
if(!params.containsKey('bowtie2_index')){
    params.bowtie2_index = null
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome R64-1-1 -profile singularity"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

// Validate some of the more idiosyncratic parameters specific to callingcards
if (params.split_fastq_by_size != null && params.split_fastq_by_part != null){
    exit 1, 'You have specified both `split_fastq_by_size` and `split_fastq_by_part`.' +
    ' Please specify only one of these parameters. The other should be null.'
}

// Check that nonsensical combinations of parameters are not set
if (params.additional_fasta && (params.bwa_index || params.bwamem2_index || params.bowtie_index || params.bowtie2_index)) {
    exit 1, 'You have specified an additional fasta file and a genome index.' +
    ' If the genome index is not equivalent to the main fasta file,' +
    ' then omit the index and allow the pipeline to create it from' +
    ' the concatenated fasta files.'
}

if (params.datatype == "mammals" && params.r1_bc_pattern == null){
    exit 1, 'You have not specified a barcode pattern for mammalian data.' +
    ' Please specify a barcode pattern using the `r1_bc_pattern` parameter.'
}

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALLINGCARDS_YEAST   } from './workflows/callingcards_yeast'
include { CALLINGCARDS_MAMMALS } from './workflows/callingcards_mammals'

//
// WORKFLOW: Run main nf-core/callingcards analysis pipeline
//
workflow NFCORE_CALLINGCARDS_MAMMALS {
    CALLINGCARDS_MAMMALS ()
}

workflow NFCORE_CALLINGCARDS_YEAST {
    CALLINGCARDS_YEAST ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    if(params.datatype == 'mammals'){
        NFCORE_CALLINGCARDS_MAMMALS ()
    } else if (params.datatype == 'yeast'){
        NFCORE_CALLINGCARDS_YEAST ()
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
