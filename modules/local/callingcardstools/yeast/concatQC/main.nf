process CONCATQC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/callingcardstools:1.2.0--pyhdfd78af_0' :
        'biocontainers/callingcardstools:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(qc_files), path(barcode_details)

    output:
    tuple val(meta), path("*.csv") , optional: false , emit: qc
    path "versions.yml"            , optional: false , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    callingcardstools yeast_combine_qc \\
        -i ${qc_files.join(' ')} \\
        -b ${barcode_details} \\
        -p ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callingcardstools: \$(echo \$(callingcardstools --version 2>&1) | sed 's/.* //')
    END_VERSIONS

    """
}
