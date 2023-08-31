#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FETCHTOOL_ASSEMBLY } from '../../../../../modules/ebi-metagenomics/fetchtool/assembly/main.nf'

workflow test_fetchtool_assembly {

    input = [
        [ id:'ERZ1022727' ], // meta map
        "ERZ1022727",
    ]

    config = file("tests/modules/ebi-metagenomics/fetchtool/configs/testing.json")

    FETCHTOOL_ASSEMBLY ( input, config )
}
