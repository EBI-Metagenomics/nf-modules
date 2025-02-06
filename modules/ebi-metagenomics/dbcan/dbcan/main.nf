process DBCAN {
    tag "$meta.id"
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
    def args = task.ext.args ?: ''
    def dbcan_version = "${params.dbcan_version}"
    def db_version = "${params.db_version}"
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
    def dbcan_version = "${params.dbcan_version}"
    def db_version = "${params.db_version}"
    """
    mkdir dbcan
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: ${dbcan_version}
        dbcan_db: ${db_version}
    END_VERSIONS
    """
}
