
include { GENERATEGAF                              } from '../../../modules/ebi-metagenomics/generategaf/main'
include { OWLTOOLS                                 } from '../../../modules/ebi-metagenomics/owltools/main'
include { MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS  } from '../../../modules/ebi-metagenomics/mgnifypipelinestoolkit/summarisegoslims/main'

workflow GOSLIM_SWF {

    take:
    ch_ips     // channel: [ val(meta), path(tsv) ]
    go_obo     // file: path(obo)
    goslim_ids // file: path(txt)
    go_banding // file: path(txt)

    main:

    ch_versions = Channel.empty()

    GENERATEGAF( ch_ips )
    ch_versions = ch_versions.mix(GENERATEGAF.out.versions.first())

    OWLTOOLS(
        GENERATEGAF.out.gaf,
        go_obo,
        goslim_ids
    )
    ch_versions = ch_versions.mix(OWLTOOLS.out.versions.first())

    // Sync channels
    ch_ips.join( OWLTOOLS.out.gaf )
        .multiMap { meta, interproscan_tsv, gaf_out ->
            interproscan: [meta, interproscan_tsv ]
            gaf: [meta, gaf_out ]
        }.set {
            summarise_ch
        }

    MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS(
        summarise_ch.interproscan,
        summarise_ch.gaf,
        go_obo,
        go_banding
    )
    ch_versions = ch_versions.mix(MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS.out.versions.first())

    emit:
    go_summary            = MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS.out.go_summary     // channel: [ val(meta), path(csv) ]
    goslim_summary        = MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS.out.goslim_summary // channel: [ val(meta), path(csv) ]
    versions              = ch_versions                                                // channel: [ versions.yml ]
}

