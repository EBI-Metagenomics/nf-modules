
process MAPSEQ2BIOM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(msq)
    tuple path(db_fasta), path(db_tax), path(db_otu), path(db_mscluster), val(db_label)

    output:
    tuple val(meta), path("${meta.id}.txt")                                 , emit: krona_input
    tuple val(meta), path("${meta.id}.tsv"), path("${meta.id}.notaxid.tsv") , emit: biom_out
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    mapseq2biom \
        $args
        --krona ${prefix}.txt \
        --no-tax-id-file ${prefix}.notaxid.tsv \
        --label $label \
        --query ${mapseq_out} \
        --otu-table ${db_otu} \
        --out-file ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq2biom: 0.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    touch ${prefix}.notaxid.tsv
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq2biom: 0.1.0
    END_VERSIONS
    """
}
