process CONCATQC {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmatkhan/callingcardstools:0.2.0' :
        'docker://cmatkhan/callingcardstools:0.2.0' }"

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
