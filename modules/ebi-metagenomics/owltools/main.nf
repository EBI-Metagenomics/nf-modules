
process OWLTOOLS {
    tag "$meta.id"
    label 'process_single'

    container "microbiome-informatics/owltools:2024-06-12"

    input:
    tuple val(meta), path(input_gaf)
    path go_obo
    path goslim_ids

    output:
    tuple val(meta), path("*.gaf"), emit: gaf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    owltools \\
        ${go_obo} \\
        --gaf ${input_gaf} \\
        $args \\
        --idfile ${goslim_ids} \\
        --write-gaf ${prefix}_goslim_annotations.gaf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        owltools: 2024-06-12
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_goslim_annotations.gaf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        owltools: 2024-06-12
    END_VERSIONS
    """
}
