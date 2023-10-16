#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTA_DOMAINANNOTATION } from '../../../../subworkflows/ebi-metagenomics//main.nf'

workflow test_fasta_domainannotation {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    FASTA_DOMAINANNOTATION ( input )
}
