process CHECKM2_DOWNLOAD_DB {

    label 'process_single'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    output:
    path("out/CheckM2_database/*.dmnd")                 , emit: checkm2_db
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    checkm2 database --download --path out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CheckM2 : \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir -p out/CheckM2_database
    touch out/CheckM2_database/uniref100.KO.1.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CheckM2 : \$(checkm2 --version)
    END_VERSIONS
    """
}
