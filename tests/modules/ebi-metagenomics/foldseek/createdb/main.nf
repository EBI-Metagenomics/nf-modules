#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDSEEK_CREATEDB } from '../../../../../modules/ebi-metagenomics/foldseek/createdb/main.nf'

workflow test_foldseek_createdb {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['proteomics']['pdb']['tim1_pdb'], checkIfExists: true) ]
    ]

    FOLDSEEK_CREATEDB ( input )
}
