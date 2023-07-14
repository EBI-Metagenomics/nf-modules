#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INFERNAL_CMSEARCH } from '../../../../../modules/ebi-metagenomics/infernal/cmsearch/main.nf'

workflow test_infernal_cmsearch {
    
    input = [
        [ id:'test', single_end:false ],                                                               // meta map
        file("tests/modules/ebi-metagenomics/infernal/cmsearch/data/example.fa", checkIfExists: true), // seqdb
    ]

    cm_model = file("tests/modules/ebi-metagenomics/infernal/cmsearch/data/tRNA5.c.cm", checkIfExists: true)

    INFERNAL_CMSEARCH ( input , cm_model )
}
