#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BMTAGGER } from '../../../../../modules/ebi-metagenomics/bmtagger/bmtagger/main.nf'
include { CREATE_DB_BMTAGGER } from '../../../../../modules/ebi-metagenomics/bmtagger/index_reference/main.nf'

workflow test_bmtagger {

    input = [
        [ id:'test_fasta', single_end:true ],                                                               // meta map
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/example.fa", checkIfExists: true), // input
    ]

    reference_fasta = file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true)

    CREATE_DB_BMTAGGER(reference_fasta)

    input_format = channel.value("fasta")
    output_directory = channel.value("bmtagger_output")

    BMTAGGER ( input, CREATE_DB_BMTAGGER.out.bitmask, CREATE_DB_BMTAGGER.out.srprism, input_format, output_directory )
}

workflow test_bmtagger_fastq {

    input = [
        [ id:'test_fasta', single_end:true ],                                                               // meta map
        [file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_1.fastq", checkIfExists: true),
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_2.fastq", checkIfExists: true)
        ], // input
    ]

    reference_fasta = file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true)

    CREATE_DB_BMTAGGER(reference_fasta)

    input_format = channel.value("fastq")
    output_directory = channel.value("bmtagger_output")

    BMTAGGER ( input, CREATE_DB_BMTAGGER.out.bitmask, CREATE_DB_BMTAGGER.out.srprism, input_format, output_directory )
}
