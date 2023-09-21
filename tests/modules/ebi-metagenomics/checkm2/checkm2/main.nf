#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM2               } from '../../../../../modules/ebi-metagenomics/checkm2/checkm2/main.nf'
include { CHECKM2_DOWNLOAD_DB   } from '../../../../../modules/ebi-metagenomics/checkm2/download_db/main.nf'

workflow test_checkm2 {
    meta = [ id:'test', single_end:false ]

    CHECKM2_DOWNLOAD_DB()

    input = [
        meta,
        file("./tests/modules/ebi-metagenomics/checkm2/checkm2/data/bins", checkIfExists: true)
    ]

    checkm_db = CHECKM2_DOWNLOAD_DB.out.checkm2_db

    CHECKM2 ( input, checkm_db)
}

workflow test_checkm2_empty_directory {
    meta = [ id:'test', single_end:false ]

    input = [
        meta,
        file("./tests/modules/ebi-metagenomics/checkm2/checkm2/data/empty_bins/*")
    ]

    checkm_db = file("./tests/modules/ebi-metagenomics/checkm2/checkm2/data/NO_FILE")

    CHECKM2 (input, checkm_db)
}
