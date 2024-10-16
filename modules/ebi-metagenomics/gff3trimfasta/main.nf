process GFF3_TRIM_FASTA {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/gawk:4.1.3--0'

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*._trimmed.gff"), optional: true, emit: gff
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """

    awk '/##FASTA/{exit}1' "$tab" > "${prefix}_trimmed.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | grep -o '[0-9]\{8\}')
    END_VERSIONS
    """

    stub:
    """
    touch "${prefix}_trimmed.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | grep -o '[0-9]\{8\}')
    END_VERSIONS
    """
}