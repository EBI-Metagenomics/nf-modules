process GFF2GBK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcbio-gff:0.7.1--pyhdfd78af_3':
        'biocontainers/bcbio-gff:0.7.1--pyhdfd78af_3' }"

    input:
    tuple val(meta), path(contigs), path(gff), path(proteins)

    output:
    tuple val(meta), path("*_from_gff.gbk"), emit: gbk
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_arg = proteins.size() > 0 ? "--proteins ${proteins}" : ""
    """
    gbk_generator.py \\
        $args \\
        --contigs ${contigs} \\
        --gff ${gff} \\
        --output_gbk ${prefix}_from_gff.gbk \\
        ${proteins_arg}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_arg = proteins.size() > 0 ? "--proteins ${proteins}" : ""
    """
    echo $args

    touch ${prefix}_from_gff.gbk
    """
}
