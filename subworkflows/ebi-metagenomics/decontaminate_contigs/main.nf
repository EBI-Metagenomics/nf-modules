include { MINIMAP2_ALIGN } from '../../../modules/ebi-metagenomics/minimap2/align/main'
include { FILTERPAF      } from '../../../modules/ebi-metagenomics/filterpaf/main'
include { SEQKIT_GREP    } from '../../../modules/ebi-metagenomics/seqkit/grep/main'

workflow DECONTAMINATE_CONTIGS {
    take:
    contigs            // [ meta, path(assembly_fasta)]
    contaminant_genome // [ meta, path(reference_fasta)]

    main:
    ch_versions = Channel.empty()

    // Mix and match the genomes and contigs
    contigs
        .join(contaminant_genome)
        .multiMap { meta, assembly, cont_genome ->
            contigs: [meta, assembly]
            genome: [meta, cont_genome]
        }
        .set { synced_contigs_contaminant_genome }

    MINIMAP2_ALIGN(
        synced_contigs_contaminant_genome.contigs,
        synced_contigs_contaminant_genome.genome,
        false,             // bam_format
        false,             // bam_index_extension
        false,             // cigar_paf_format
        false              // cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    FILTERPAF(
        MINIMAP2_ALIGN.out.paf
    )
    ch_versions = ch_versions.mix(FILTERPAF.out.versions.first())

    contigs
        .join(FILTERPAF.out.mapped_contigs_txt)
        .multiMap { meta, assembly, mapped_contigs_txt ->
            contigs: [meta, assembly]
            pattern: mapped_contigs_txt
        }
        .set { seqkit_grep_input_ch }

    SEQKIT_GREP(
        seqkit_grep_input_ch.contigs,
        seqkit_grep_input_ch.pattern,
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    cleaned_contigs = SEQKIT_GREP.out.filter
    versions        = ch_versions
}
