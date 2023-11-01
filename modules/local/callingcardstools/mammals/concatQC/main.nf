process CONCATQC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/callingcardstools:1.2.0--pyhdfd78af_0' :
        'biocontainers/callingcardstools:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(qbed_files)
    tuple val(meta), path(barcode_qc_files)

    output:
    tuple val(meta), path("*.qbed")            , optional: false , emit: qbed
    tuple val(meta), path("*_barcode_qc.tsv")  , optional: false , emit: barcode_qc
    tuple val(meta), path("*_srt_count.tsv")   , optional: false , emit: srt_count
    tuple val(meta), path("*_aln_summary.tsv") , optional: false , emit: aln_summary
    path "versions.yml"                        , optional: false , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.prefix ?: ""

    """
    callingcardstools mammals_combine_qc \\
        -q ${qbed_files.join(' ')} \\
        -b ${barcode_qc_files.join(' ')} \\
        -f ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callingcardstools: \$(echo \$(callingcardstools --version 2>&1) | sed 's/.* //')
    END_VERSIONS

    """
}
