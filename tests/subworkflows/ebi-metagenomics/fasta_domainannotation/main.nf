#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/ebi-metagenomics/diamond/makedb/main.nf'
include { FASTA_DOMAINANNOTATION } from '../../../../subworkflows/ebi-metagenomics//main.nf'

workflow test_fasta_domainannotation {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    blast_db = BLAST_MAKEBLASTDB ( fasta ).out.db

    eggnog_db = file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/eggnog.db", checkIfExists: true)
    eggnog_data_dir = eggnog_db.parent
    diamond_db = DIAMOND_MAKEDB ( fasta ).out.db
    ch_eggnog = [ eggnog_db, eggnog_data_dir, diamond_db ]

    FASTA_DOMAINANNOTATION ( [ [id:'test'], fasta ], blast_db, ch_eggnog )
}
