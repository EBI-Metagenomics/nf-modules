#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FRAGGENESCAN } from '../../../../modules/ebi-metagenomics/fraggenescan/main.nf'

workflow test_fraggenescan {

    input = [
        [ id:'test', single_end:false ],                                                               // meta map
        file("tests/modules/ebi-metagenomics/fraggenescan/data/example.fa.gz", checkIfExists: true), // seqdb
    ]


    FRAGGENESCAN ( input )
}
