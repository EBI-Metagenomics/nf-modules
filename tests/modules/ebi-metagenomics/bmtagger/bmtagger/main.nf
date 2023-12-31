#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BMTAGGER_BMTAGGER } from '../../../../../modules/ebi-metagenomics/bmtagger/bmtagger/main.nf'
include { BMTAGGER_INDEXREFERENCE } from '../../../../../modules/ebi-metagenomics/bmtagger/indexreference/main.nf'

workflow test_bmtagger_single_end {

    // TODO the following tests are disabled, bmtagger takes way too long to run in the CI server.

    input = [
        [ id:'test_fasta', single_end:true ], // meta map
        [
            file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_1.fastq", checkIfExists: true)
        ], // input
    ]
    input_ref = [
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true),
    ]

    BMTAGGER_INDEXREFERENCE( input_ref )

    input_format = channel.value("fastq")
    output_directory = channel.value("bmtagger_output")

    BMTAGGER_BMTAGGER (
        input,
        BMTAGGER_INDEXREFERENCE.out.bitmask,
        BMTAGGER_INDEXREFERENCE.out.srprism,
        input_format,
        output_directory
    )
}

workflow test_bmtagger_paired_end {

    // TODO the following tests are disabled, bmtagger takes way too long to run in the CI server.

    input = [
        [ id:'test_fasta', paired_end:true ], // meta map
        [
            file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_1.fastq", checkIfExists: true),
            file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/ex_2.fastq", checkIfExists: true)
        ], // input
    ]
    input_ref = [
        file("tests/modules/ebi-metagenomics/bmtagger/bmtagger/data/reference.fasta", checkIfExists: true),
    ]

    BMTAGGER_INDEXREFERENCE( input_ref )

    input_format = channel.value("fastq")
    output_directory = channel.value("bmtagger_output")

    BMTAGGER_BMTAGGER (
        input,
        BMTAGGER_INDEXREFERENCE.out.bitmask,
        BMTAGGER_INDEXREFERENCE.out.srprism,
        input_format,
        output_directory
    )
}
