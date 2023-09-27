#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SEQ } from '../../../../../modules/ebi-metagenomics/seqtk/seq/main.nf'

workflow test_seqtk_seq {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/seqtk/seq/data/SRR21814853.merged.fastq.gz', checkIfExists: true)
    ]

    SEQTK_SEQ ( input )
}
