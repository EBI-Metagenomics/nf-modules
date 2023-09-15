#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTP } from '../../../../modules/ebi-metagenomics/fastp/main.nf'

workflow test_fastp {
    
    input_1 = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
        ]

    input_2 = true

    input_3 = true

    FASTP (
        input_1,
        input_2,
        input_3
     )
}
