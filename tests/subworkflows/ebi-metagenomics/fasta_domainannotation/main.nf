#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/ebi-metagenomics/diamond/makedb/main.nf'
include { FASTA_DOMAINANNOTATION } from '../../../../subworkflows/ebi-metagenomics/fasta_domainannotation/main.nf'

workflow test_fasta_domainannotation {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    input = Channel.of( [ [id:'test'], fasta ] )

    blast_fasta = Channel.value( fasta )

    eggnog_db = Channel.value( file("tests/modules/ebi-metagenomics/eggnogmapper/data/fixtures/eggnog.db", checkIfExists: true) )
    eggnog_data_dir = eggnog_db.parent
    diamond_db = DIAMOND_MAKEDB ( fasta ).db

    FASTA_DOMAINANNOTATION ( input, blast_fasta, eggnog_db, eggnog_data_dir, diamond_db )
}
