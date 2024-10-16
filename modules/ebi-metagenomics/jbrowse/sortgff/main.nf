process SORT_GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/coreutils:8.25--0'

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("${meta.id}_sorted.gff"), optional: true, emit: gff
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    (grep "^#" $tab; grep -v "^#" $tab | sort -t"\$(printf '\\t')" -k1,1 -k4,4n) > ${prefix}_sorted.gff;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}_sorted.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
