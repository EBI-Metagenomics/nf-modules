#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EASEL_ESLSFETCH } from '../../../../../modules/ebi-metagenomics/easel/eslsfetch/main.nf'

workflow test_easel_eslsfetch {
    
    input_fasta = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/easel/eslsfetch/data/test.fasta.gz', checkIfExists: true)
    ]

    input_cmsearchdeoverlap = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/easel/eslsfetch/data/test.tblout.deoverlapped', checkIfExists: true)
    ]

    EASEL_ESLSFETCH ( input_fasta, input_cmsearchdeoverlap )
}
