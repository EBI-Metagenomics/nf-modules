process BGCSMAPPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/622d8944750bc95bb56b4c3ed5c2b827e677c14073d48a5231e0f2bec0718add/data':
        'biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(gff), path(sanntis), path(gecco), path(antismash)

    output:
    tuple val(meta), path("*_bgcs.gff"), optional: true, emit: gff
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sanntis_param = sanntis ? "--sanntis_gff ${sanntis}" : ""
    def gecco_param = gecco ? "--gecco_gff ${gecco}" : ""
    def antismash_param = antismash ? "--antismash_gff ${antismash}" : ""
    """
    bgc_mapper.py \\
        ${sanntis_param} \\
        ${gecco_param} \\
        ${antismash_param} \\
        --base_gff ${gff} \\
        --output_gff ${prefix}_bgcs.gff
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_bgcs.gff
    """
}
