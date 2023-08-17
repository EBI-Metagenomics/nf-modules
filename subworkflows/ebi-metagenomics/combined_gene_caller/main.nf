
include { PRODIGAL                 } from '../../../modules/nf-core/prodigal/main'
include { FRAGGENESCAN             } from '../../../modules/ebi-metagenomics/fraggenescan/main'
include { COMBINEDGENECALLER_MERGE } from '../../../modules/ebi-metagenomics/combinedgenecaller/merge/main'

workflow  COMBINED_GENE_CALLER {

    take:
    ch_assembly // channel: [ val(meta), [ fasta ] ]
    ch_mask_file // channel: [ .out | .txt ]

    main:

    ch_versions = Channel.empty()

    PRODIGAL ( ch_assembly )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions.first())

    FRAGGENESCAN ( ch_assembly )
    ch_versions = ch_versions.mix(FRAGGENESCAN.out.versions.first())

    ch_annotations = PRODIGAL.out.

    COMBINEDGENECALLER_MERGE (  )

        tuple val(meta), path(prodigal_out), path(prodigal_ffn), path(prodigal_faa), path(fgs_out), path(fgs_ffn), path(fgs_faa)
    path(mask_file)

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

