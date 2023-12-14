#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RRNA_EXTRACTION } from '../../../../subworkflows/ebi-metagenomics/rrna_extraction/main.nf'

workflow test_rrna_extraction {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/subworkflows/ebi-metagenomics/rrna_extraction/data/test.fasta.gz', checkIfExists: true)
    ]

    input_ch = Channel.of( input )

    RRNA_EXTRACTION ( input_ch )
}
