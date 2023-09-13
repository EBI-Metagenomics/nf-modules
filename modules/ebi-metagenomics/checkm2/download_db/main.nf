process CHECKM2_DOWNLOAD_DB {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("*/*.dmnd"), emit: checkm2_db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    checkm2 database --download

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(checkm2 --version ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p out
    touch out/fake.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(checkm2 --version ))
    END_VERSIONS
    """
}
