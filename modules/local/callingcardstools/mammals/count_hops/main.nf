process COUNT_HOPS {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/callingcardstools:1.2.0--pyhdfd78af_0' :
        'biocontainers/callingcardstools:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(barcode_details)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*_qbed.pkl")        , optional: true , emit: qbed_pkl
    tuple val(meta), path("*_barcode_qc.pkl")  , optional: true , emit: barcode_qc_pkl
    tuple val(meta), path("*.qbed")            , optional: true , emit: qbed
    tuple val(meta), path("*_barcode_qc.tsv")  , optional: true , emit: barcode_qc
    tuple val(meta), path("*_srt_count.tsv")   , optional: true , emit: srt_count
    tuple val(meta), path("*_aln_summary.tsv") , optional: true , emit: aln_summary
    tuple val(meta), path("*passing.bam")      , optional: true , emit: passing_bam
    tuple val(meta), path('*failing.bam')      , optional: true , emit: failing_bam
    path "versions.yml"                        , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/callingcards/bin/
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: ""
    """

    callingcardstools process_mammals_bam \\
    -i $bam \\
    -b $barcode_details \\
    -g $fasta \\
    -f $prefix \\
    -s $suffix \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callingcardstools: \$(callingcardstools --version | sed 's/callingcardstools //g')
    END_VERSIONS
    """
}
