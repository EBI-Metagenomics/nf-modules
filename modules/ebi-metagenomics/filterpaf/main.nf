process FILTERPAF {
    tag "$meta.id"

    input:
    tuple val(meta), path(paf_file)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: unmapped_contigs_txt

    script:
    """
    awk '((\$4 - \$3) / \$2) > ${params.min_qcov} && \$12 > ${params.min_mapq}' ${paf_file} | cut -f 1 > ${meta.id}.txt
    """
}