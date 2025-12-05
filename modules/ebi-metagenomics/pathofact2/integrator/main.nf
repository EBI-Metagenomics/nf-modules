process PATHOFACT2_INTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(gff), path(toxins), path(virulence)

    output:
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tox_param = toxins ? "--toxins ${toxins}" : ""
    def vf_param = virulence ? "--virulence ${virulence}" : ""
    """
    pathofact2_integrator.py \\
        -c ${gff} \\
        ${tox_param} \\
        ${vf_param} \\
        -o ${prefix}_pathofact2_integrated.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pathofact2_integrated.gff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
