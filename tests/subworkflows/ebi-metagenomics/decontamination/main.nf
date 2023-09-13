#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DECONTAMINATION_WITH_BMTAGGER } from '../../../../subworkflows/ebi-metagenomics/decontamination/main.nf'

workflow test_decontamination {

    input = [
        [ id:'test_fasta', single_end:true ],                                                               // meta map
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/example.fa", checkIfExists: true), // input
    ]

    reference_fasta = file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true)

    CREATE_DB_BMTAGGER(reference_fasta)

    DECONTAMINATION_WITH_BMTAGGER ( input, Channel.empty() )
}
