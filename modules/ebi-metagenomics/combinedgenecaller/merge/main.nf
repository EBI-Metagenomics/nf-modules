
process COMBINEDGENECALLER_MERGE {
    tag "$meta.id"
    label 'process_single'

    container 'microbiome-informatics/combined-gene-caller:v1.0.4'

    input:
    tuple val(meta), path(prodigal_out), path(prodigal_ffn), path(prodigal_faa), path(fgs_out), path(fgs_ffn), path(fgs_faa)
    path(mask_file)

    output:
    tuple val(meta), path("*.faa"), emit: faa
    tuple val(meta), path("*.ffn"), emit: ffn
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    combined_gene_caller \
        -n $prefix \
        --prodigal-out $prodigal_out \
        --prodigal-ffn $prodigal_ffn \
        --prodigal-faa $prodigal_faa \
        --fgs-out $fgs_out \
        --fgs-ffn $fgs_ffn \
        --fgs-faa $fgs_faa \
        --mask $mask_file \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combined_gene_caller: \$(echo \$(combined_gene_caller --version))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa
    touch ${prefix}.ffn

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combined_gene_caller: \$(echo \$(combined_gene_caller --version))
    END_VERSIONS
    """
}
