process EUKGENERENAME {
    tag "${meta.id}"
    label 'process_low'

    container  "quay.io/biocontainers/biopython:1.76--2"

    input:
    val meta
    path(braker_aa)
    path(braker_ffn)

    output:
    path( "${meta.id}_braker_unique.faa"   ), emit: renamed_braker_aa
    path( "${meta.id}_braker_unique.ffn"   ), emit: renamed_braker_ffn


    script:
    """
    gene_rename.py \\
        --fasta $braker_aa \\
        --output ${meta.id} 
    """
}

// rename braker protein IDs suffix with |braker_only
