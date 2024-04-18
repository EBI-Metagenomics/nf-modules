
include { MAPSEQ             } from '../../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2BIOM        } from '../../../modules/ebi-metagenomics/mapseq2biom/main'
include { KRONA_KTIMPORTTEXT } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'

workflow MAPSEQ_OTU_KRONA {

    take:
    ch_fasta    // channel: [ val(meta), [ bam ] ]
    ch_dbs      // channel: [ path(fasta), path(tax), path(otu), path(mscluster), val(label) ]

    main:

    ch_versions = Channel.empty()

    db_fasta = ch_dbs[0]
    db_tax = ch_dbs[1]
    db_otu = ch_dbs[2]
    db_mscluster = ch_dbs[3]
    db_label = ch_dbs[4]

    MAPSEQ(
        ch_fasta,
        tuple(db_fasta, db_tax, db_mscluster)
    )
    ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())

    MAPSEQ2BIOM(
        MAPSEQ.out.mseq,
        tuple(db_otu, db_label)
    )
    ch_versions = ch_versions.mix(MAPSEQ2BIOM.out.versions.first())

    KRONA_KTIMPORTTEXT(
        MAPSEQ2BIOM.out.krona_input
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
    mseq                  = MAPSEQ.out.mseq                   // channel: [ val(meta), [ mseq ] ]
    krona_input           = MAPSEQ2BIOM.out.krona_input       // channel: [ val(meta), [ txt ] ]
    biom_out              = MAPSEQ2BIOM.out.biom_out          // channel: [ val(meta), [ tsv ] ]
    biom_notaxid_out      = MAPSEQ2BIOM.out.biom_notaxid_out  // channel: [ val(meta), [ tsv ] ]
    html                  = KRONA_KTIMPORTTEXT.out.html       // channel: [ val(meta), [ bam ] ]
    versions              = ch_versions                       // channel: [ versions.yml ]
}
