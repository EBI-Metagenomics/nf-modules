#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAPSEQ } from '../../../../modules/ebi-metagenomics/mapseq/main.nf'

workflow test_mapseq_silva_ssu {
    
    subunit_input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/mapseq/data/subunit_silva_SSU.fasta', checkIfExists: true)
    ]

    db_fasta = file('tests/modules/ebi-metagenomics/mapseq/data/silva_ssu/silva_ssu-20200130_SSU_trimmed.fasta', checkIfExists: true)
    db_tax = file('tests/modules/ebi-metagenomics/mapseq/data/silva_ssu/silva_ssu-20200130_SSU_filtered2.txt', checkIfExists: true)
    db_mscluster = file('tests/modules/ebi-metagenomics/mapseq/data/silva_ssu/silva_ssu-20200130_SSU_trimmed.fasta.mscluster', checkIfExists: true)

    db_tuple = tuple(db_fasta, db_tax, db_mscluster)

    MAPSEQ ( subunit_input,  db_tuple )
}

workflow test_mapseq_silva_lsu {
    
    subunit_input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/mapseq/data/subunit_silva_LSU.fasta', checkIfExists: true)
    ]

    db_fasta = file('tests/modules/ebi-metagenomics/mapseq/data/silva_lsu/silva_lsu-20200130_LSU_trimmed.fasta', checkIfExists: true)
    db_tax = file('tests/modules/ebi-metagenomics/mapseq/data/silva_lsu/silva_lsu-20200130_LSU_filtered2.txt', checkIfExists: true)
    db_mscluster = file('tests/modules/ebi-metagenomics/mapseq/data/silva_lsu/silva_lsu-20200130_LSU_trimmed.fasta.mscluster', checkIfExists: true)

    db_tuple = tuple(db_fasta, db_tax, db_mscluster)

    MAPSEQ ( subunit_input,  db_tuple )
}

workflow test_mapseq_unite {
    
    subunit_input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/mapseq/data/subunit_unite.fasta', checkIfExists: true)
    ]

    db_fasta = file('tests/modules/ebi-metagenomics/mapseq/data/unite/UNITE-20200214_unite_trimmed.fasta', checkIfExists: true)
    db_tax = file('tests/modules/ebi-metagenomics/mapseq/data/unite/UNITE-20200214_UNITE-tax.txt', checkIfExists: true)
    db_mscluster = file('tests/modules/ebi-metagenomics/mapseq/data/unite/UNITE-20200214_unite_trimmed.fasta.mscluster', checkIfExists: true)

    db_tuple = tuple(db_fasta, db_tax, db_mscluster)

    MAPSEQ ( subunit_input,  db_tuple )
}

workflow test_mapseq_itsonedb {
    
    subunit_input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/ebi-metagenomics/mapseq/data/subunit_itsonedb.fasta', checkIfExists: true)
    ]

    db_fasta = file('tests/modules/ebi-metagenomics/mapseq/data/itsonedb/ITSoneDB-20200214_itsonedb_trimmed.fasta', checkIfExists: true)
    db_tax = file('tests/modules/ebi-metagenomics/mapseq/data/itsonedb/ITSoneDB-20200214_ITSonedb-tax.txt', checkIfExists: true)
    db_mscluster = file('tests/modules/ebi-metagenomics/mapseq/data/itsonedb/ITSoneDB-20200214_itsonedb_trimmed.fasta.mscluster', checkIfExists: true)

    db_tuple = tuple(db_fasta, db_tax, db_mscluster)

    MAPSEQ ( subunit_input,  db_tuple )
}