#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EGGNOGMAPPER } from '../../../../modules/ebi-metagenomics/eggnogmapper/main.nf'

workflow test_eggnogmapper {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/test_queries.fa", checkIfExists: true)
    ]
    eggnog_data_dir = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures",                        checkIfExists: true)
    eggnog_diamond_db = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/eggnog_proteins.dmnd", checkIfExists: true)
    eggnog_db = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/eggnog.db",                    checkIfExists: true)

    EGGNOGMAPPER ( input, eggnog_data_dir, eggnog_diamond_db, eggnog_db )
}
