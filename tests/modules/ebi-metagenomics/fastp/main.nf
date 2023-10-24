#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTP } from '../../../../modules/ebi-metagenomics/fastp/main.nf'

workflow test_fastp {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
        ]

    FASTP (
        input,
        params.save_trimmed_fail,
        true // save merged reads
     )
}
