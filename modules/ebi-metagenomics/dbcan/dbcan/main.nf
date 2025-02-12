process DBCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/dbcan:4.1.4--pyhdfd78af_0'
        : 'biocontainers/dbcan:4.1.4--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(faa), path(gff)
    tuple path(dbcan_db, stageAs: "dbcan_db"), val(db_version)

    output:
    tuple val(meta), path("dbcan/", type: "dir"), emit: dbcan_output
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    run_dbcan \\
        ${args} \\
        --dia_cpu ${task.cpus} \\
        --hmm_cpu ${task.cpus} \\
        --tf_cpu ${task.cpus} \\
        --dbcan_thread ${task.cpus} \\
        --db_dir \$(pwd)/dbcan_db \\
        --out_dir dbcan \\
        --cgc_substrate \\
        --cluster ${gff} \\
        ${faa} \\
        protein

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 4.1.4
        dbcan_db: "${db_version}"
    END_VERSIONS
    """

    stub:
    """
    mkdir dbcan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 4.1.4
        dbcan_db: ${db_version}
    END_VERSIONS
    """
}
