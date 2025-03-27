
include { GENERATEGAF                              } from '../../../modules/ebi-metagenomics/generategaf/main'
include { OWLTOOLS                                 } from '../../../modules/ebi-metagenomics/owltools/main'
include { MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS  } from '../../../modules/ebi-metagenomics/mgnifypipelinestoolkit/summarisegoslims/main'

workflow GOSLIM_SWF {

    take:
    ch_ips     // channel: [ val(meta), path(tsv) ]
    go_obo     // channel: path(obo)
    goslim_ids // channel: path(txt)
    go_banding // channel: path(txt)

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

    MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS(
        ch_ips,
        OWLTOOLS.out.gaf,
        go_obo,
        go_banding
    )
    ch_versions = ch_versions.mix(MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS.out.versions.first())

    emit:
    go_summary            = MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS.out.go_summary     // channel: [ val(meta), path(csv) ]
    goslim_summary        = MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS.out.goslim_summary // channel: [ val(meta), path(csv) ]
    versions              = ch_versions                                                // channel: [ versions.yml ]
}

