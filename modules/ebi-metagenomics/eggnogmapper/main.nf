process EGGNOGMAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::eggnog-mapper=2.1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(eggnog_data_dir)
    path(eggnog_diamond_db)
    path(eggnog_db)

    output:
    tuple val(meta), path("*.emapper.hits"), emit: csv
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.12'
    """
    emapper.py \\
        --cpu ${task.cpus} \\
        -i ${fasta} \\
        --data_dir ${eggnog_data_dir} \\
        -m diamond \\
        --dmnd_db ${eggnog_diamond_db} \\
        --database ${eggnog_db} \\
        --output ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.12'
    """
    touch ${prefix}.emapper.hits

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: $VERSION
    END_VERSIONS
    """
}
