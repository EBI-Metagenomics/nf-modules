process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*._sorted.gff"), optional:true, emit: gff
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """

    (grep "^#" $tab; grep -v "^#" $tab | sort -t"`printf '\t'`" -k1,1 -k4,4n)  > ${prefix}.sorted.gff;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version| awk '{print $NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.sorted.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version| awk '{print $NF}')
    END_VERSIONS
    """
}
