
process SUMMARISEGOSLIMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.2.0--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ips)
    tuple val(meta), path(gaf)
    path go_obo
    path go_banding

    output:
    tuple val(meta), path("*_summary.csv"), emit: go_summary
    tuple val(meta), path("*_slim.csv")   , emit: goslim_summary
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    summarise_goslims \\
        -go ${go_obo} \\
        -gb ${go_banding} \\
        -i ${ips} \\
        -gaf ${gaf} \\
        -o ${prefix}_summary

    mv ${prefix}_summary ${prefix}_summary.csv
    mv ${prefix}_summary_slim ${prefix}_summary_slim.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.csv
    touch ${prefix}_summary_slim.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
