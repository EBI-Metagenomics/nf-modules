#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { READ_QC } from '../../../../subworkflows/ebi-metagenomics/read_qc/main.nf'

workflow test_read_qc {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
        ]

    READ_QC ( input )
}
