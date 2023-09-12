#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EGGNOGMAPPER } from '../../../../modules/ebi-metagenomics/eggnogmapper/main.nf'

workflow test_eggnogmapper {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    EGGNOGMAPPER ( input )
}
