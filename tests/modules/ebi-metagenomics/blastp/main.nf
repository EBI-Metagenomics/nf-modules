#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLASTP } from '../../../../modules/ebi-metagenomics/blastp/main.nf'

workflow test_blastp {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("tests/modules/ebi-metagenomics/blastp/data/test.fa", checkIfExists: true)
    ]
    db = file("tests/modules/ebi-metagenomics/blastp/data/uniprot_sprot_subset_DB", checkIfExists: true)

    BLASTP ( input, db )
}
