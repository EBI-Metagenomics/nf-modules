process PATHOFACT2_EXTRACTFASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84':
        'biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(fasta), path(blastp_out), path(pathofact2_tox), path(pathofact2_vf)

    output:
    tuple val(meta), path("*_pathofact2.fasta"), optional: true, emit: fasta
    tuple val(meta), path("*_support.tsv")     , optional: true, emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pathofact_fasta_extractor.py \\
        -f ${fasta} \\
        -b ${blastp_out} \\
        -t ${pathofact2_tox} \\
        -v ${pathofact2_vf} \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}_pathofact2.fasta
    touch ${prefix}_support.tsv
    """
}
