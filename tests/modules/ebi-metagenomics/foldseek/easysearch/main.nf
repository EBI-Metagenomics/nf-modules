#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDSEEK_CREATEDB   } from '../../../../../modules/ebi-metagenomics/foldseek/createdb/main.nf'
include { FOLDSEEK_EASYSEARCH } from '../../../../../modules/ebi-metagenomics/foldseek/easysearch/main.nf'

workflow test_foldseek_easysearch {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("tests/modules/ebi-metagenomics/foldseek/createdb/data", checkIfExists: true)
    ]

    input2 = [
        [ id:'test2', single_end:false ], // meta map
        file("tests/modules/ebi-metagenomics/foldseek/createdb/data/1tim.pdb", checkIfExists: true)
    ]

    FOLDSEEK_CREATEDB ( input )
    FOLDSEEK_EASYSEARCH ( input2, FOLDSEEK_CREATEDB.out.db )
}
