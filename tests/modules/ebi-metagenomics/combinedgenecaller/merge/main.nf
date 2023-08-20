#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COMBINEDGENECALLER_MERGE } from '../../../../../modules/ebi-metagenomics/combinedgenecaller/merge/main.nf'

workflow test_combinedgenecaller_merge {

    input = [
        [ id:'test', single_end:false ],                                                           // meta map
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.out.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.ffn.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.faa.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.out.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.ffn.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.faa.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/cmsearch.all.tblout.deoverlapped.gz", checkIfExists: true)
    ]

    COMBINEDGENECALLER_MERGE ( input )
}


workflow test_combinedgenecaller_merge_no_mask {

    input = [
        [ id:'test', single_end:false ],                                                           // meta map
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.out.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.ffn.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_prodigal.faa.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.out.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.ffn.gz", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/combinedgenecaller/merge/data/input_fgs.faa.gz", checkIfExists: true),
        null // mask file
    ]

    COMBINEDGENECALLER_MERGE ( input )
}
