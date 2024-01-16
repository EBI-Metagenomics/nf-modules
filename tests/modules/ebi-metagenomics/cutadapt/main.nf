#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUTADAPT } from '../../../../modules/ebi-metagenomics/cutadapt/main.nf'

workflow test_cutadapt_pe_double_trim {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4822665_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4822665_2.fastq.gz', checkIfExists: true) ],
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/fwd_primer.fasta', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/rev_primer.fasta', checkIfExists: true) ]
        ]

    CUTADAPT ( input )
}

workflow test_cutadapt_pe_fwd_trim {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4822665_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4822665_2.fastq.gz', checkIfExists: true) ],
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/fwd_primer.fasta', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/empty.txt') ]
        ]

    CUTADAPT ( input )
}

workflow test_cutadapt_pe_rev_trim {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4822665_1.fastq.gz', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4822665_2.fastq.gz', checkIfExists: true) ],
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/empty.txt'),
          file('tests/modules/ebi-metagenomics/cutadapt/data/rev_primer.fasta', checkIfExists: true) ]
        ]

    CUTADAPT ( input )
}

workflow test_cutadapt_se_double_trim {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4334351.fastq.gz', checkIfExists: true)],
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/fwd_primer.fasta', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/rev_primer.fasta', checkIfExists: true) ]
        ]

    CUTADAPT ( input )
}

workflow test_cutadapt_se_fwd_trim {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4334351.fastq.gz', checkIfExists: true)],
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/fwd_primer.fasta', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/empty.txt', checkIfExists: true) ]
        ]

    CUTADAPT ( input )
}

workflow test_cutadapt_se_rev_trim {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/ERR4334351.fastq.gz', checkIfExists: true)],
        [ file('tests/modules/ebi-metagenomics/cutadapt/data/empty.txt', checkIfExists: true),
          file('tests/modules/ebi-metagenomics/cutadapt/data/rev_primer.fasta', checkIfExists: true) ]
        ]

    CUTADAPT ( input )
}