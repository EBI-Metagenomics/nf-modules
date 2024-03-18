#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ESMFOLD } from '../../../../modules/ebi-metagenomics/esmfold/main.nf'

workflow test_esmfold {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/esmfold/data/one_protein.fasta", checkIfExists: true) // params.test_data['sarscov2']['genome']['proteome_fasta']
    ]

    ESMFOLD ( input, "cpu" )
}

workflow test_esmfold_gpu {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/esmfold/data/one_protein.fasta", checkIfExists: true)
    ]

    ESMFOLD ( input, "gpu" )
}
