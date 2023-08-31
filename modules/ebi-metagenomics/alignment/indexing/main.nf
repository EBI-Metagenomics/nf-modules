process INDEX_FASTA {

    tag "${meta.id} index ${fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${fasta.baseName}*.?*"), emit: fasta_with_index

    script:
    """
    echo "index ref genome"
    bwa-mem2 index ${fasta}
    """
}
