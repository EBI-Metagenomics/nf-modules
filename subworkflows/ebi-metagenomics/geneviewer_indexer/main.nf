
nextflow.enable.dsl = 2

include { INDEX_FASTA } from './index_fasta/index_fasta'

Channel.of(params.input).set { inputChannel }
outputChannel = params.output_dir

workflow {
    INDEX_FASTA(inputChannel,outputChannel)
}

