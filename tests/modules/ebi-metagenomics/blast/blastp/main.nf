#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_BLASTP } from '../../../../../modules/ebi-metagenomics/blast/blastp/main.nf'

workflow test_blastp {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("tests/modules/ebi-metagenomics/blast/blastp/data/test.fa", checkIfExists: true)
    ]
    db = file("tests/modules/ebi-metagenomics/blast/blastp/data/uniprot_sprot_subset_DB", checkIfExists: true)

    BLAST_BLASTP ( input, db )
}
