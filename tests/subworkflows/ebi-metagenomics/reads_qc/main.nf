#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { READS_QC } from '../../../../subworkflows/ebi-metagenomics/reads_qc/main.nf'

workflow test_reads_qc {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/fastp/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
        ]

    READS_QC ( input )
}
