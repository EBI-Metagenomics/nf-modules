
include { MAPSEQ             } from '../../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2BIOM        } from '../../../modules/ebi-metagenomics/mapseq2biom/main'
include { KRONA_KTIMPORTTEXT } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'

workflow MAPSEQ_OTU_KRONA {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fasta    // channel: [ val(meta), [ bam ] ]
    ch_dbs      // channel: [ path(fasta), path(tax), path(otu), path(mscluster), val(label) ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    db_fasta = ch_dbs[0]
    db_tax = ch_dbs[1]
    db_otu = ch_dbs[2]
    db_mscluster = ch_dbs[3]
    db_label = ch_dbs[4]
    

    MAPSEQ(
        ch_fasta,
        tuple(db_fasta, db_tax, db_mscluster)
    )


    // SAMTOOLS_SORT ( ch_bam )
    // ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    mseq      = MAPSEQ.out.mseq          // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    // versions = ch_versions                     // channel: [ versions.yml ]
}

