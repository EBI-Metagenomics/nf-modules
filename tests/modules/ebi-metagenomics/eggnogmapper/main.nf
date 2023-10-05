#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/ebi-metagenomics/diamond/makedb/main.nf'
include { EGGNOGMAPPER } from '../../../../modules/ebi-metagenomics/eggnogmapper/main.nf'

workflow test_eggnogmapper {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    eggnog_data_dir = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures",                        checkIfExists: true)
    // eggnog_diamond_db = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/eggnog_proteins.dmnd", checkIfExists: true)
    eggnog_db = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/eggnog.db",                    checkIfExists: true)

    DIAMOND_MAKEDB ( fasta )
    EGGNOGMAPPER ( [ [id:'test'], fasta ], eggnog_data_dir, DIAMOND_MAKEDB.out.db, eggnog_db )
}
