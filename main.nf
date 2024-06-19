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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALLINGCARDS  } from './workflows/callingcards'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_callingcards_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_callingcards_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_callingcards_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
params.gtf = getGenomeAttribute('gtf')

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
        params.bwa_index = getGenomeAttribute(params, 'bwa')
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
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CALLINGCARDS {

    take:
    reads
    barcode_details

    main:

    if(params.fasta){
        ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).collect()
            .map{ it -> [[id:it[0].getSimpleName()], it[0]]}
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

    //
    // WORKFLOW: Run pipeline
    //
    CALLINGCARDS (
        reads,
        barcode_details,
        ch_fasta,
        ch_gtf,
        ch_regions_mask,
        additional_fasta,
        rseqc_modules
    )

    emit:
    multiqc_report = CALLINGCARDS.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CALLINGCARDS (
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.barcode_details
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CALLINGCARDS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
