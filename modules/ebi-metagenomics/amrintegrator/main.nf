process AMRINTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.2.11--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.2.11--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(deeparg), path(rgi), path(amrfp), path(gff)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    amrintegrator \\
        --deeparg_out ${deeparg} \\
        --rgi_out ${rgi} \\
        --amrfp_out ${amrfp} \\
        --cds_gff ${gff} \\
        --output ${prefix}.gff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
