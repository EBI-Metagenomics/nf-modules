process DBCAN {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:4.1.4--pyhdfd78af_0' :
        'biocontainers/dbcan' }"

    input:
    tuple val(meta), path(faa), path(gff)
    tuple path(dbcan_db, stageAs: "dbcan_db"), val(db_version)

    output:
    tuple path("dbcan/", type: "dir"), val("${params.db_version}"), emit: dbcan_output
    path "versions.yml"                                      , emit: versions

    script:
    """
    run_dbcan \\
        --dia_cpu ${task.cpus} \\
        --hmm_cpu ${task.cpus} \\
        --tf_cpu ${task.cpus} \\
        --dbcan_thread ${task.cpus} \\
        --db_dir dbcan_db \\
        --out_dir dbcan \\
        --cgc_substrate \\
        --cluster ${gff} \\
        ${faa} \\
        protein

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: "${params.dbcan_version}"
        dbcan_db: "${params.db_version}"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = "${params.dbcan_version}"
    def DB_VERSION = "${params.db_version}"
    """
    mkdir dbcan
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: ${VERSION}
        dbcan_db: ${DB_VERSION}
    END_VERSIONS
    """
}
