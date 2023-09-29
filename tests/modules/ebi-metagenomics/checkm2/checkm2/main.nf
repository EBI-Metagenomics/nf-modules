#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM2               } from '../../../../../modules/ebi-metagenomics/checkm2/checkm2/main.nf'
include { CHECKM2_DOWNLOAD_DB   } from '../../../../../modules/ebi-metagenomics/checkm2/download_db/main.nf'

workflow test_checkm2 {

    CHECKM2_DOWNLOAD_DB()
    checkm_db = CHECKM2_DOWNLOAD_DB.out.checkm2_db

    meta = [ id:'test', single_end:false ]
    input = [
        meta,
        [file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    ]
    CHECKM2 ( input, checkm_db )
}

workflow test_checkm2_empty_directory {

    meta = [ id:'test', single_end:false ]
    // TODO: change file to params.test_data when test DB comes
    input = [
        meta,
        file("./tests/modules/ebi-metagenomics/checkm2/checkm2/data/*.fa")
    ]

    // TODO: change to dmnd DB when it appears
    checkm_db = file("./tests/modules/ebi-metagenomics/checkm2/checkm2/data/NO_FILE")

    CHECKM2 ( input, checkm_db )
}
