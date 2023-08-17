
process COMBINEDGENECALLER_MERGE {
    tag "$meta.id"
    label 'process_single'

    container 'microbiome-informatics/combined-gene-caller:v1.0.4'

    input:
    tuple val(meta), path(prodigal_out), path(prodigal_ffn), path(prodigal_faa), path(fgs_out), path(fgs_ffn), path(fgs_faa)
    path(mask_file)

    output:
    tuple val(meta), path("*.faa.gz"), emit: faa
    tuple val(meta), path("*.ffn.gz"), emit: ffn
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf $prodigal_out > prodigal_out_uncompressed
    gzip -cdf $prodigal_ffn > prodigal_ffn_uncompressed
    gzip -cdf $prodigal_faa > prodigal_faa_uncompressed
    gzip -cdf $fgs_out > fgs_out_uncompressed
    gzip -cdf $fgs_ffn > fgs_ffn_uncompressed
    gzip -cdf $fgs_faa > fgs_faa_uncompressed
    gzip -cdf $mask_file > mask_file_uncompressed

    combined_gene_caller \\
        -n $prefix \\
        --prodigal-out prodigal_out_uncompressed \\
        --prodigal-ffn prodigal_ffn_uncompressed \\
        --prodigal-faa prodigal_faa_uncompressed \\
        --fgs-out fgs_out_uncompressed \\
        --fgs-ffn fgs_ffn_uncompressed \\
        --fgs-faa fgs_faa_uncompressed \\
        --mask mask_file_uncompressed \\
        $args

    gzip -n $prefix.*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combined_gene_caller: \$(combined_gene_caller --version)
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
