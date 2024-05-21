process CONCATFASTA {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/coreutils:8.25--1' :
    'biocontainers/coreutils:8.25--1' }"

    input:
    tuple val(meta), path(main_fasta)
    path(additional_fasta)

    output:
    tuple val(meta), path("*concat.fasta"), emit: fasta
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: ""
    def suffix      = task.ext.suffix   ?: ""
    def filename    = prefix + suffix + "concat.fasta"
    def VERSION     = "8.25"    // WARN: Version information not provided by
                                // tool on CLI. Please update this string
                                // when bumping container versions.

    """
    cat ${main_fasta} ${additional_fasta} > ${filename}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
