#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTMHMM } from '../../../../modules/ebi-metagenomics/deeptmhmm/main.nf'

workflow test_deeptmhmm {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    DEEPTMHMM ( fasta )
}
