#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/ebi-metagenomics/diamond/makedb/main.nf'
include { FASTA_DOMAINANNOTATION } from '../../../../subworkflows/ebi-metagenomics/fasta_domainannotation/main.nf'

workflow test_fasta_domainannotation {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    input = Channel.of( [ [id:'test'], fasta ] )
    blast_fasta = Channel.value( fasta )
    blast_mode = "diamond" // "diamond" or "blast"

    FASTA_DOMAINANNOTATION ( input, blast_fasta, blast_mode )
}
