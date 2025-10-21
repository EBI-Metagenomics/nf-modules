process EUKGENEOVERLAP {
    tag "${meta.id}"
    label 'process_medium'

    container  "quay.io/biocontainers/biopython:1.76--2"

    input:
    val meta
    path(braker_gff)
    path(metaeuk_gff)
    path(braker_aa)
    path(metaeuk_aa)
    val threshold

    output:
    path( "${meta.id}_overlap.tsv"         ), emit: overlap_tsv
    path( "${meta.id}_braker_unique.tsv"   ), emit: braker_unique
    path( "${meta.id}_metaeuk_unique.tsv"  ), emit: metaeuk_unique
    path( "${meta.id}_overlap.faa"         ), emit: overlap_faa
    path( "${meta.id}_braker_unique.faa"   ), emit: braker
    path( "${meta.id}_metaeuk_unique.faa"  ), emit: metaeuk


    script:
    def threshold_args     = threshold ? "--threshold ${threshold}": '--threshold 0.9'

    """
    gene_overlap.py \\
        --braker_gff $braker_gff \\
        --metaeuk_gff $metaeuk_gff \\
        --braker_aa $braker_aa \\
        --metaeuk_aa $metaeuk_aa \\
        --output ${meta.id} \\
        $threshold_args
    """
}
