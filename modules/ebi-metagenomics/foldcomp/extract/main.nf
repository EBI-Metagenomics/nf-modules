process FOLDCOMP_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldcomp:0.0.5--h43eeafb_2':
        'biocontainers/foldcomp:0.0.5--h43eeafb_2' }"

    input:
    tuple val(meta), path(fcz)
    val(mode)

    output:
    tuple val(meta), path("{${fcz}_plddt,*.plddt}"), optional: true, emit: plddt
    tuple val(meta), path("{${fcz}_fasta,*.fasta}"), optional: true, emit: fasta
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode_arg = mode == "fasta" ? "--fasta" : "--plddt"
    """
    foldcomp \\
        extract \\
        -t ${task.cpus} \\
        ${mode_arg} \\
        ${fcz} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldcomp: v0.0.5
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plddt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldcomp: v0.0.5
    END_VERSIONS
    """
}
