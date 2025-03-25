process MGNIFYPIPELINESTOOLKIT_KRONATXTFROMCATCLASSIFICATION {
    tag "${meta.id}"
    label 'process_single'

    container 'microbiome-informatics/mgnify-pipelines-toolkit:1.0.3'

    input:
    tuple val(meta), path(cat_output)
    tuple val(meta2), path(taxonomy)

    output:
    tuple val(meta), path("${meta.id}.krona.txt"), emit: krona_txt
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    krona_txt_from_cat_classification \\
        --input ${cat_output} \\
        --output ${meta.id}.krona.txt \\
        --names_dmp ${taxonomy}/names.dmp \\
        --nodes_dmp ${taxonomy}/nodes.dmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.krona.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
