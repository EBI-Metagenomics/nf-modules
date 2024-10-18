process TAXONKIT_REFORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxonkit:0.17.0--h9ee0642_1':
        'biocontainers/taxonkit:0.17.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(tsv)
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: reformat_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    taxonkit \\
        reformat \\
        $args \\
        --threads $task.cpus \\
        --data-dir $taxdb \\
        --out-file ${prefix}.tsv \\
        $tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$( taxonkit version | sed 's/.* v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$( taxonkit version | sed 's/.* v//' )
    END_VERSIONS
    """
}
