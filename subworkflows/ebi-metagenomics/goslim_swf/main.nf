
include { GENERATEGAF           } from '../../../modules/ebi-metagenomics/generategaf/main'
include { OWLTOOLS              } from '../../../modules/ebi-metagenomics/owltools/main'
include { SUMMARISEGOSLIMS      } from '../../../modules/ebi-metagenomics/summarisegoslims/main'

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

    OWLTOOLS ( 
        GENERATEGAF.out.gaf,
        go_obo,
        goslim_ids
    )
    ch_versions = ch_versions.mix(OWLTOOLS.out.versions.first())

    SUMMARISEGOSLIMS (
        ch_ips,
        OWLTOOLS.out.gaf,
        go_obo,
        go_banding
    )
    ch_versions = ch_versions.mix(SUMMARISEGOSLIMS.out.versions.first())

    emit:
    go_summary            = SUMMARISEGOSLIMS.out.go_summary     // channel: [ val(meta), path(csv) ]
    goslim_summary        = SUMMARISEGOSLIMS.out.goslim_summary // channel: [ val(meta), path(csv) ]
    versions              = ch_versions                         // channel: [ versions.yml ]
}

