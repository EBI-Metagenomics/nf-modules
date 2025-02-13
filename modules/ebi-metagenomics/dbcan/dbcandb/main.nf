process DBCANDB {

    tag "dbCAN $params.db_version"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val(db_ftp_link)
    val(db_version)

    output:
    tuple path("dbcan_db/", type: "dir"), val("${params.db_version}"), emit: dbcan_db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget ${db_ftp_link}

    tar -xvzf dbcan_"${db_version}".tar.gz

    # delete the next line when using the right db version!
    mv 4.1.3-V12 dbcan_db

    rm dbcan_"${db_version}".tar.gz

    touch dbcan_db/dbcan_"${db_version}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan_db: "${db_version}"
    END_VERSIONS
    """

    stub:
    """
    mkdir dbcan_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan_db: "${db_version}"
    END_VERSIONS
    """
}
