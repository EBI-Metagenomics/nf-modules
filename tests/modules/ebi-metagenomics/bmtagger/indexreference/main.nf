#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BMTAGGER_INDEXREFERENCE } from '../../../../../modules/ebi-metagenomics/bmtagger/indexreference/main.nf'

workflow test_bmtagger_indexreference {

    input = [
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true),
    ]

    BMTAGGER_INDEXREFERENCE ( input )
}
