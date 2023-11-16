#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDCOMP_COMPRESS   } from '../../../../../modules/ebi-metagenomics/foldcomp/compress/main.nf'
include { FOLDCOMP_DECOMPRESS } from '../../../../../modules/ebi-metagenomics/foldcomp/decompress/main.nf'

workflow test_foldcomp_decompress_dir {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/foldcomp/compress/data/", checkIfExists: true)
    ]
    FOLDCOMP_COMPRESS ( input )
    FOLDCOMP_DECOMPRESS ( FOLDCOMP_COMPRESS.out.fcz )
}

workflow test_foldcomp_decompress_single {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/foldcomp/compress/data/4590051235.pdb", checkIfExists: true)
    ]
    FOLDCOMP_COMPRESS ( input )
    FOLDCOMP_DECOMPRESS ( FOLDCOMP_COMPRESS.out.fcz )
}
