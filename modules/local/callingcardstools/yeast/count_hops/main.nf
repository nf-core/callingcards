process COUNT_HOPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/callingcardstools:1.2.0--pyhdfd78af_0' :
        'biocontainers/callingcardstools:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(barcode_details)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.qbed")                                                 , optional: false , emit: qbed
    tuple val(meta), path("*_summary.tsv")                                          , optional: false , emit: qc_summary
    tuple val(meta), path("*_passing_tagged.bam"), path("*_passing_tagged.bam.bai") , optional: false , emit: passing_bam
    tuple val(meta), path("*_failing_tagged.bam"), path("*_failing_tagged.bam.bai") , optional: false , emit: failing_bam
    path "versions.yml"                                                             , optional: false , emit: versions
    tuple val(meta), path("*_aln_info.tsv")                                         , optional: true  , emit: aln_info

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    callingcardstools process_yeast_bam \\
        -i ${bam} \\
        -g ${fasta} \\
        -j ${barcode_details} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callingcardstools: \$(echo \$(callingcardstools --version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
