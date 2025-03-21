
process CONVERTCMSCANTOCMSEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.0.2--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cmscan_tblout)

    output:
    tuple val(meta), path("*cmsearch.tbl"), emit: cmsearch_tblout
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    convert-cmscan-to-cmsearch-tblout.py --input $cmscan_tblout --output ${meta.id}.cmsearch.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cmsearch.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
