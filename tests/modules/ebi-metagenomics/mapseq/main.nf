#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAPSEQ } from '../../../../modules/ebi-metagenomics/mapseq/main.nf'

workflow test_mapseq {
    
    subunit_input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/mapseq/data/SRR17062740_SSU_trimmed.fasta', checkIfExists: true)
    ]

    db_fasta = file('tests/modules/ebi-metagenomics/mapseq/data/silva/silva_ssu-20200130_SSU_trimmed.fasta', checkIfExists: true)
    db_tax = file('tests/modules/ebi-metagenomics/mapseq/data/silva/silva_ssu-20200130_SSU_filtered2.txt', checkIfExists: true)
    db_mscluster = file('tests/modules/ebi-metagenomics/mapseq/data/silva/silva_ssu-20200130_SSU_trimmed.fasta.mscluster', checkIfExists: true)

    db_tuple = tuple(db_fasta, db_tax, db_mscluster)

    MAPSEQ ( subunit_input,  db_tuple )
}
