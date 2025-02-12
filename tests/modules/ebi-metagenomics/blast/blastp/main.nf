#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../../modules/ebi-metagenomics/blast/makeblastdb/main.nf'
include { BLAST_BLASTP } from '../../../../../modules/ebi-metagenomics/blast/blastp/main.nf'

workflow test_blast_blastp {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    out_ext = '' // empty test case to check default
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}

workflow test_blast_blastp_xml {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    out_ext = 5
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}

workflow test_blast_blastp_tsv {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    out_ext = 102
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}

workflow test_blast_blastp_csv {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    out_ext = 48
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}
