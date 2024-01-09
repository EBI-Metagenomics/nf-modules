process COLABFOLD_COLABFOLDBATCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/colabfold:1.5.5--pyh7cba7a3_0':
        'biocontainers/colabfold:1.5.5--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}"), emit: out_dir
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    colabfold_batch \\
        $args \\
        ${fasta} \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: 1.5.2
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: 1.5.2
    END_VERSIONS
    """
}
