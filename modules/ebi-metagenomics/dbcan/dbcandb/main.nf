process DBCANDB {
    tag "$meta.id"
    tag "dbCAN $params.db_version"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    output:
    tuple val(meta), path("dbcan_db/", type: "dir"), val("${params.db_version}"), emit: dbcan_db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


script:
    def args = task.ext.args ?: ''
    def db_ftp_link = "${params.db_ftp_link}"
    def db_version = "${params.db_version}"
    """
    wget ${db_ftp_link}

    tar -xvzf dbcan_"${db_version}".tar.gz

    # delete the next line when using the right db version!
    mv ${db_version} dbcan_db

    rm dbcan_"${db_version}".tar.gz

    touch dbcan_db/dbcan_"${db_version}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan_db: "${db_version}"
    END_VERSIONS
    """
    
stub:
    def args = task.ext.args ?: ''
    def	db_version = "${db_version}"
    """
    mkdir dbcan_db
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan_db: "${db_version}"
    END_VERSIONS
    """
}
