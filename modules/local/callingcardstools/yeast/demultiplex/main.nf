process DEMULTIPLEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/callingcardstools:1.2.0--pyhdfd78af_0' :
        'biocontainers/callingcardstools:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(barcode_details)

    output:
    tuple val(meta), path("*.fq")     , optional: false , emit: reads
    path "versions.yml"               , optional: false , emit: versions
    tuple val(meta), path("*.pickle") , optional: true  , emit: pickle
    tuple val(meta), path("*.csv")    , optional: true  , emit: qc
    tuple val(meta), path("*.txt")    , optional: true  , emit: verbose_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    callingcardstools split_fastq \\
        -r1 ${reads[0]} \\
        -r2 ${reads[1]} \\
        -b ${barcode_details} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        callingcardstools: \$(echo \$(callingcardstools --version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
