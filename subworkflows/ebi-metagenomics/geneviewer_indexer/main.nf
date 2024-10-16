nextflow.enable.dsl = 2

include { INDEX_FASTA } from './index_fasta/index_fasta'
include { INDEX_GFF } from './index_gff/index_gff'

// Set the input and output channels using parameters from config
Channel.of(params.fasta_input).set { fastaInputChannel }
Channel.of(params.gff_input).set { gffInputChannel }

fastaOutputChannel = params.fasta_output_dir
gffOutputChannel = params.gff_output_dir

workflow {
    // Run both subworkflows
    INDEX_FASTA(fastaInputChannel, fastaOutputChannel)
    INDEX_GFF(gffInputChannel, gffOutputChannel)
}
