process SANNTIS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/sanntis:0.9.3.5--pyhdfd78af_1' :
        'biocontainers/sanntis:0.9.3.5--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(interproscan_output), path(gbk)

    output:
    tuple val(meta), path("*_sanntis.gff"), emit: gff
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -v "/protein_id=" ${gbk} > ${prefix}_prepped.gbk
    sanntis \
        --ip-file ${interproscan_output} \
        --outfile ${prefix}_sanntis.gff \
        --cpu ${task.cpus} \
        ${args} \
        ${prefix}_prepped.gbk


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sanntis.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
    END_VERSIONS
    """
}
