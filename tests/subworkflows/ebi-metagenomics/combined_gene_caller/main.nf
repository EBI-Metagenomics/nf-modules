#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COMBINED_GENE_CALLER } from '../../../../subworkflows/ebi-metagenomics/combined_gene_caller/main.nf'

workflow test_combined_gene_caller {

    input = [
        [ id:'test', single_end:false ],
        file("tests/subworkflows/ebi-metagenomics/combined_gene_caller/data/ERZ19591644_FASTA_subsample.fasta.gz", checkIfExists: true)
    ]

    COMBINED_GENE_CALLER ( input, Channel.empty() )
}
