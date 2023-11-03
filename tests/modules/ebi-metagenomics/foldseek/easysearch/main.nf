#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDSEEK_EASYSEARCH } from '../../../../../modules/ebi-metagenomics/foldseek/easysearch/main.nf'

workflow test_foldseek_easysearch {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    FOLDSEEK_EASYSEARCH ( input )
}
