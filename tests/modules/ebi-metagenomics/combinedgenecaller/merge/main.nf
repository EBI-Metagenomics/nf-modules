#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COMBINEDGENECALLER_MERGE } from '../../../../../modules/ebi-metagenomics/combinedgenecaller/merge/main.nf'

workflow test_combinedgenecaller_merge {

    input = [
        [ id:'test', single_end:false ],                                                           // meta map
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.out", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.ffn", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.faa", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.out", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.ffn", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.faa", checkIfExists: true)
    ]
    mask_file = file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/cmsearch.all.tblout.deoverlapped", checkIfExists: true)

    COMBINEDGENECALLER_MERGE ( input, mask_file )
}
