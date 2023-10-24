#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { READS_QC } from '../../../../subworkflows/ebi-metagenomics/reads_qc/main.nf'

workflow test_reads_qc_pe {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
    ]

    READS_QC ( input, false )
}

workflow test_reads_qc_se {
    
    input = [
        [ id:'test', single_end:true ], // meta map
          file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR9674626.fastq.gz', checkIfExists: true),
    ]

    READS_QC ( input, false  )
}

workflow test_reads_qc_pe_and_se {
    
    input_pe = [
        [ id:'test_pe', single_end:false ], // meta map
        [ file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
    ]

    input_se = [
        [ id:'test_se', single_end:true ], // meta map
          file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR9674626.fastq.gz', checkIfExists: true)
    ]

    input_pe_and_se = Channel.from( input_pe, input_se )

    READS_QC ( input_pe_and_se, false  )
}

workflow test_save_merged {
    
    input_pe = [
        [ id:'test_pe', single_end:false ], // meta map
        [ file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR21814853_1.fastq.gz', checkIfExists: true),
          file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR21814853_2.fastq.gz', checkIfExists: true) ]
    ]

    input_se = [
        [ id:'test_se', single_end:true ], // meta map
          file('tests/subworkflows/ebi-metagenomics/reads_qc/data/SRR9674626.fastq.gz', checkIfExists: true)
    ]

    input_pe_and_se = Channel.from( input_pe, input_se )

    READS_QC ( input_pe_and_se, true )
}
