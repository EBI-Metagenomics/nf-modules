#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {  } from '../../../../subworkflows/ebi-metagenomics//main.nf'

workflow test_combined_gene_caller {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

     ( input )
}
