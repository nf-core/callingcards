process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads), path(barcode_details)

    output:
    path '*R1.fq'                        , emit: r1
    path barcode_details                 , emit: barcode_details
    tuple val(meta), path '*.tsv'        , emit: id_to_bc_map
    path "versions.yml"                  , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/callingcards/bin/
    """
    parse_fastq.py \\
        -r1 $reads[0] \\
        -r2 $reads[1] \\
        -b $reads[3] \\
        -o $component_id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
