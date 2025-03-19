process KRONA_TXT_FROM_CAT_CLASSIFICATION {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.1--9856f872fdeac74e':
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), path(cat_output)

    output:
    tuple val(meta), path("${meta.id}.krona.txt"), emit: krona_txt
    path "versions.yml"                          , emit: versions

    script:
    """
    krona_txt_from_CAT_classification.py \\
        -i ${cat_output} \\
        -o ${meta.id}.krona.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}