include { BLAST_MAKEBLASTDB } from '../../../modules/ebi-metagenomics/blast/makeblastdb/main'
include { BLAST_BLASTP      } from '../../../modules/ebi-metagenomics/blast/blastp/main'
include { INTERPROSCAN      } from '../../../modules/ebi-metagenomics/interproscan/main'
include { EGGNOGMAPPER      } from '../../../modules/ebi-metagenomics/eggnogmapper/main'

workflow FASTA_DOMAINANNOTATION {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

