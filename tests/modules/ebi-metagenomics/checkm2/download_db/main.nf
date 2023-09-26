#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM2_DOWNLOAD_DB   } from '../../../../../modules/ebi-metagenomics/checkm2/download_db/main.nf'

workflow test_checkm2_download_db {

    CHECKM2_DOWNLOAD_DB()

}
