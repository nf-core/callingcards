/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME                                   } from "../subworkflows/local/prepare_genome.nf"
include { PREPARE_READS as YEAST_PREPARE_READS             } from "../subworkflows/local/yeast/prepare_reads.nf"
include { PREPARE_READS as MAMMALS_PREPARE_READS           } from "../subworkflows/local/mammals/prepare_reads.nf"
include { ALIGN                                            } from "../subworkflows/local/align.nf"
include { PROCESS_ALIGNMENTS as YEAST_PROCESS_ALIGNMENTS   } from "../subworkflows/local/yeast/process_alignments.nf"
include { PROCESS_ALIGNMENTS as MAMMALS_PROCESS_ALIGNMENTS } from "../subworkflows/local/mammals/process_alignments.nf"
include { MULTIQC                                          } from "../modules/nf-core/multiqc"
include { paramsSummaryMap                                 } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                           } from '../subworkflows/local/utils_nfcore_callingcards_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CALLINGCARDS {

    take:
    ch_reads
    ch_barcode_details
    ch_fasta
    ch_gtf
    ch_regions_mask
    additional_fasta
    rseqc_modules

    main:

    ch_versions          = Channel.empty()
    ch_multiqc_files     = Channel.empty()
    ch_fasta_index       = Channel.empty()
    ch_bam_index         = Channel.empty()
    ch_prepared_reads    = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flatstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()


    // SUBWORKFLOW_1: Index the genome in various ways
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

    // SUBWORKFLOW_2: Prepare the reads (trimming, removing callingcards specific barcodes)
    if(params.datatype == 'yeast'){
        YEAST_PREPARE_READS (
            ch_reads,
            ch_barcode_details
        )

        ch_versions = ch_versions.mix(YEAST_PREPARE_READS.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PREPARE_READS.out.raw_fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PREPARE_READS.out.trimmomatic_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PREPARE_READS.out.demux_fastqc_zip.collect{it[1]}.ifEmpty([]))

        ch_prepared_reads = ch_prepared_reads.mix(YEAST_PREPARE_READS.out.reads)

    } else if (params.datatype == 'mammals'){
        MAMMALS_PREPARE_READS (
            ch_reads
        )
        ch_versions = ch_versions.mix(MAMMALS_PREPARE_READS.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PREPARE_READS.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PREPARE_READS.out.trimmomatic_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PREPARE_READS.out.umi_log.collect{it[1]}.ifEmpty([]))

        ch_prepared_reads = ch_prepared_reads.mix(MAMMALS_PREPARE_READS.out.reads)
    }

    ALIGN (
        ch_prepared_reads,
        PREPARE_GENOME.out.bwamem2_index,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.bowtie_index
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    if (params.datatype == 'yeast'){
        // create a channel from the ch_barcode details
        // of structure [ [id: sample_id], path(barcode_details) ]
        // this will be used to join the barcode_details to the alignments
        ch_barcode_details
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
        // SUBWORKFLOW_3: process the alignnment into a qbed file
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
        YEAST_PROCESS_ALIGNMENTS (
            ch_aln_with_details,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.genome_bed,
            rseqc_modules,
            ch_gtf
        )
        ch_versions = ch_versions.mix(YEAST_PROCESS_ALIGNMENTS.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.samtools_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.samtools_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.samtools_idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.picard_qc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.rseqc_bamstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.rseqc_inferexperiment.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.rseqc_innerdistance.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.rseqc_readdistribution.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.rseqc_readduplication.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(YEAST_PROCESS_ALIGNMENTS.out.rseqc_tin.collect{it[1]}.ifEmpty([]))


    } else if (params.datatype == 'mammals'){
        // join the alignment output with the ch_barcode_details
        // to create a channel with structure:
        // [val(meta), path(bam), path(bai), path(barcode_details)]
        ALIGN.out.bam
            .map{meta, bam, bai -> [meta.id, meta, bam, bai] }
            .combine(ch_barcode_details
                        .map{meta, barcode_details ->
                                [meta.id, barcode_details]},
                    by: 0)
            .map{id, meta, bam, bai, barcode_details ->
                    [meta, bam, bai, barcode_details] }
            .set{ ch_aln_with_details }

        //
        // SUBWORKFLOW_4: process the alignnment into a qbed file
        //
        MAMMALS_PROCESS_ALIGNMENTS (
            ch_aln_with_details,
            ch_fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.genome_bed,
            rseqc_modules,
            ch_gtf
        )
        ch_versions = ch_versions.mix(MAMMALS_PROCESS_ALIGNMENTS.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.samtools_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.samtools_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.samtools_idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.picard_qc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.rseqc_bamstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.rseqc_inferexperiment.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.rseqc_innerdistance.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.rseqc_readdistribution.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.rseqc_readduplication.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MAMMALS_PROCESS_ALIGNMENTS.out.rseqc_tin.collect{it[1]}.ifEmpty([]))

    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
