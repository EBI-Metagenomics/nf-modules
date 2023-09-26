#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BMTAGGER_INDEXREFERENCE } from '../../../../../modules/ebi-metagenomics/bmtagger/indexreference/main.nf'

workflow test_bmtagger_indexreference {

    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BMTAGGER_INDEXREFERENCE ( input )
}
