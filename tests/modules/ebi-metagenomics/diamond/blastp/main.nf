#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../../modules/ebi-metagenomics/diamond/makedb/main.nf'
include { DIAMOND_BLASTP } from '../../../../../modules/ebi-metagenomics/diamond/blastp/main.nf'

workflow test_diamond_blastp {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    out_ext = 'txt'
    blast_columns = 'qseqid qlen'

    DIAMOND_MAKEDB ( [ [id:'test'], fasta ] )
    DIAMOND_BLASTP ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db.map { meta, dmnd -> dmnd }, out_ext, blast_columns )
}

workflow test_diamond_blastp_daa {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    out_ext = 'daa'
    blast_columns = []

    DIAMOND_MAKEDB ( [ [id:'test'], fasta ] )
    DIAMOND_BLASTP ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db.map { meta, dmnd -> dmnd }, out_ext, blast_columns )
}
