#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDCOMP_COMPRESS } from '../../../../../modules/ebi-metagenomics/foldcomp/compress/main.nf'
include { FOLDCOMP_EXTRACT  } from '../../../../../modules/ebi-metagenomics/foldcomp/extract/main.nf'

workflow test_foldcomp_extract_plddt_single {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/foldcomp/compress/data/4590785447.pdb", checkIfExists: true)
    ]
    FOLDCOMP_COMPRESS ( input )

    FOLDCOMP_EXTRACT ( FOLDCOMP_COMPRESS.out.fcz, "plddt" )
}

workflow test_foldcomp_extract_fasta_directory {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/foldcomp/compress/data", checkIfExists: true)
    ]
    FOLDCOMP_COMPRESS ( input )

    FOLDCOMP_EXTRACT ( FOLDCOMP_COMPRESS.out.fcz, "fasta" )
}
