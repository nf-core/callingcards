process PARSEBAM {
    tag "$meta.id"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmatkhan/pycallingcards:latest' :
        'docker://cmatkhan/pycallingcards:latest' }"

    input:
    tuple val(meta), path(bam), path(bai), path(barcode_details)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*qbed")       , emit: qbed
    tuple val(meta), path("*tsv")         , emit: bc_qc
    tuple val(meta), path("*passing.bam"), emit: passing_bam
    tuple val(meta), path('*failing.bam'), emit: failing_bam
    path "versions.yml"                  , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/callingcards/bin/
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    pycallingcards parse_bam \
    -i $bam \
    -b $barcode_details \
    -g $fasta \
    -l info ${args} 2> parse_bam.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycallingcards: \$(pycallingcards --version | sed 's/callingcardstools //g')
    END_VERSIONS
    """
}
