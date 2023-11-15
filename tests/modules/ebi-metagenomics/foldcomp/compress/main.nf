#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDCOMP_COMPRESS } from '../../../../../modules/ebi-metagenomics/foldcomp/compress/main.nf'

workflow test_foldcomp_compress {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/foldcomp/compress/data/", checkIfExists: true)
    ]

    FOLDCOMP_COMPRESS ( input )
}

workflow test_foldcomp_compress_single {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/foldcomp/compress/data/4590051235.pdb", checkIfExists: true)
    ]

    FOLDCOMP_COMPRESS ( input )
}
