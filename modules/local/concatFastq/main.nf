process CONCATFASTQ {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/coreutils:8.25--1' :
    'biocontainers/coreutils:8.25--1' }"

    input:
        tuple val(meta), file(reads)

    output:
        tuple val(meta), path("*concat.fastq*"), emit: reads
        path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def suffix      = task.ext.suffix   ?: ""
    def filename    = prefix + "_" + suffix + "concat.fastq"
    def VERSION     = "8.25"    // WARN: Version information not provided by
                                // tool on CLI. Please update this string
                                // when bumping container versions.
    if (params.gzip_concatenated_fastq) {
    """
    cat ${reads.join(' ')} > ${filename} && \\
    gzip ${filename}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
    } else {
    """
    cat ${reads.join(' ')} > ${filename}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
    }
}
