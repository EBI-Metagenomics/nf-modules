#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CMSEARCHTBLOUTDEOVERLAP } from '../../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main.nf'

workflow test_cmsearchtbloutdeoverlap {

    sample = [
        [ id:'test', single_end: false ], // meta map
        file("tests/modules/ebi-metagenomics/cmsearchtbloutdeoverlap/data/1.cmscan.clan.tblout", checkIfExists: true),
    ]
    claninfo = file("tests/modules/ebi-metagenomics/cmsearchtbloutdeoverlap/data/ribo.claninfo", checkIfExists: true)

    CMSEARCHTBLOUTDEOVERLAP ( sample, claninfo )
}

workflow test_cmsearchtbloutdeoverlap_decompress {

    sample = [
        [ id:'test', single_end: false ], // meta map
        file("tests/modules/ebi-metagenomics/cmsearchtbloutdeoverlap/data/1.cmscan.clan.tblout.gz", checkIfExists: true),
    ]
    claninfo = file("tests/modules/ebi-metagenomics/cmsearchtbloutdeoverlap/data/ribo.claninfo", checkIfExists: true)

    CMSEARCHTBLOUTDEOVERLAP ( sample, claninfo )
}