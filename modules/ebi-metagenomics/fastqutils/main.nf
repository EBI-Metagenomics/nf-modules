
process FASTQUTILS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq_utils:0.25.2--h50ea8bc_1':
        'biocontainers/fastq_utils:0.25.2--h50ea8bc_1' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.txt"), emit: sanity
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq_info \\
        ${fastq} \\
        2>&1 >/dev/null \\
        | tail -1 \\
        > ${prefix}_sanity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqutils: \$(fastq_info 2>&1 >/dev/null | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sanity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqutils: \$(fastq_info 2>&1 >/dev/null | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
