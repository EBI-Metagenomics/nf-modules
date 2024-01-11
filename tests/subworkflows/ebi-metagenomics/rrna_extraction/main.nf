#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RRNA_EXTRACTION } from '../../../../subworkflows/ebi-metagenomics/rrna_extraction/main.nf'

workflow test_rrna_extraction {
    
    fa = [
        [ id:'test', single_end:false ], // meta map
        file('tests/subworkflows/ebi-metagenomics/rrna_extraction/data/test.fasta.gz', checkIfExists: true)
    ]
    input_fa = Channel.of( fa )
    
    RRNA_EXTRACTION (
        input_fa,
        file(params.rfam),
        file(params.rfam_clan)
    )
}
