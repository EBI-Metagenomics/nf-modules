#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BMTAGGER } from '../../../../../modules/ebi-metagenomics/bmtagger/bmtagger/main.nf'
include { BMTAGGER_INDEX_REFERENCE } from '../../../../../modules/ebi-metagenomics/bmtagger/index_reference/main.nf'

workflow test_bmtagger {

    input = [
        [ id:'test_fasta', single_end:true ],                                                               // meta map
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/example.fa", checkIfExists: true), // input
    ]

    reference_fasta = file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true)

    BMTAGGER_INDEX_REFERENCE(reference_fasta)

    input_format = channel.value("fasta")
    output_directory = channel.value("bmtagger_output")

    BMTAGGER ( input, BMTAGGER_INDEX_REFERENCE.out.bitmask, BMTAGGER_INDEX_REFERENCE.out.srprism, input_format, output_directory )
}

workflow test_bmtagger_fastq {

    input = [
        [ id:'test_fasta', single_end:true ],                                                               // meta map
        [file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_1.fastq", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_2.fastq", checkIfExists: true)
        ], // input
    ]

    reference_fasta = file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true)

    BMTAGGER_INDEX_REFERENCE(reference_fasta)

    input_format = channel.value("fastq")
    output_directory = channel.value("bmtagger_output")

    BMTAGGER ( input, BMTAGGER_INDEX_REFERENCE.out.bitmask, BMTAGGER_INDEX_REFERENCE.out.srprism, input_format, output_directory )
}
