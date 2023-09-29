process BMTAGGER_BMTAGGER {

    tag "$meta.id"

    label 'process_low'

    debug true

    conda "bioconda::bmtagger=3.101"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bmtagger:3.101--h470a237_4':
        'biocontainers/bmtagger:3.101--h470a237_4' }"

    input:
    tuple val(meta), path(input)
    path reference_bitmask
    path reference_srprism
    val input_format
    val output_file

    output:
    tuple val(meta), path("${output_file}"), emit: output_host
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (input_format == 'fasta') {
        args += ' -q0 '
    } else if (input_format == 'fastq') {
        args += ' -q1 '
    } else {
        error "Invalid format: ${input_format}"
    }

    if (meta.single_end) {
        args += " -1 ${input} "
    } else {
        if (input[0].name.contains("_1")) {
            args += " -1 ${input[0]} -2 ${input[1]} "
        } else {
            args += " -1 ${input[1]} -2 ${input[0]} "
        }
    }
    def srprism = reference_srprism[0].baseName
    """
    bmtagger.sh -b ${reference_bitmask} \
                -x ${srprism} \
                -o ${output_file} \
                ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bmtagger: \$(bmtagger.sh -hV 2>1 | grep version | grep -Eo '[0-9.]+')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bmtagger: \$(bmtagger.sh -hV 2>1 | grep version | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
