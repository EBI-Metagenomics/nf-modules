include { MINIMAP2_ALIGN   } from '../../../modules/ebi-metagenomics/minimap2/align/main'
include { FILTER_PAF       } from '../../../modules/ebi-metagenomics/filterpaf/main'
include { SEQKIT_GREP      } from '../../../modules/ebi-metagenomics/seqkit/grep/main'

workflow DECONTAMINATE_CONTIGS {
    take:
    contigs           // [ meta, path(assembly_fasta)]
    host_genome       // [ meta, path(reference_fasta)]

    main:
    ch_versions = Channel.empty()

    MINIMAP2_ALIGN(
        contigs,
        host_genome,
        false,        // bam_format
        false,        // bam_index_extension
        false,        // cigar_paf_format
        false         // cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    FILTER_PAF(MINIMAP2_ALIGN.out.paf)
    ch_versions = ch_versions.mix(FILTER_PAF.out.versions.first())

    SEQKIT_GREP(contigs.join(FILTER_PAF.out.mapped_contigs_txt))
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    cleaned_contigs = SEQKIT_GREP.out.filter
    versions        = ch_versions

}
