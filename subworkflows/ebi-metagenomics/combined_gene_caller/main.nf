
include { PRODIGAL                 } from '../../../modules/ebi-metagenomics/prodigal/main'
include { FRAGGENESCANRS           } from '../../../modules/ebi-metagenomics/fraggenescanrs/main'
include { COMBINEDGENECALLER_MERGE } from '../../../modules/ebi-metagenomics/combinedgenecaller/merge/main'

workflow COMBINED_GENE_CALLER {

    take:
    ch_assembly  // channel: [ val(meta), [ fasta ] ]
    ch_mask_file // channel: [ val(meta), [.out | .txt] ]

    main:

    ch_versions = Channel.empty()

    PRODIGAL ( ch_assembly, channel.value("sco") )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    FRAGGENESCANRS ( ch_assembly, params.fraggenescanrs_train_model )
    ch_versions = ch_versions.mix(FRAGGENESCANRS.out.versions)

    ch_mask_file = ch_mask_file ?: Channel.empty()

    ch_annotations = PRODIGAL.out.gene_annotations.join(
        PRODIGAL.out.nucleotide_fasta,
    ).join(
        PRODIGAL.out.amino_acid_fasta,
    ).join(
        FRAGGENESCANRS.out.gene_annotations,
    ).join(
        FRAGGENESCANRS.out.nucleotide_fasta,
    ).join(
        FRAGGENESCANRS.out.amino_acid_fasta,
    ).join(
        ch_mask_file, remainder: true
    )

    ch_merge = ch_annotations.map { meta, prodigal_sco, prodigal_nt_fasta, prodigal_aa_fasta, frag_gene_annot, frag_nt_fasta, frag_aa_fasta, mask ->
        return [
            meta,
            prodigal_sco,
            prodigal_nt_fasta,
            prodigal_aa_fasta,
            frag_gene_annot,
            frag_nt_fasta,
            frag_aa_fasta,
            mask ?: []
        ]
    }

    COMBINEDGENECALLER_MERGE( ch_merge )

    ch_versions = ch_versions.mix(COMBINEDGENECALLER_MERGE.out.versions)

    emit:
    faa      = COMBINEDGENECALLER_MERGE.out.faa  // channel: [ val(meta), [ faa ] ]
    ffn      = COMBINEDGENECALLER_MERGE.out.ffn  // channel: [ val(meta), [ ffn ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}
