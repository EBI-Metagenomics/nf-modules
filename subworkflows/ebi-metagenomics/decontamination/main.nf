include { ALIGNMENT             } from '../../../modules/ebi-metagenomics/alignment/mapping/main'
include { INDEX_FASTA           } from '../../../modules/ebi-metagenomics/alignment/indexing/main'
include { SAMTOOLS_FASTQ        } from '../../../modules/nf-core/samtools/fastq/main'
include { BMTAGGER              } from '../../../modules/ebi-metagenomics/bmtagger/main'


workflow DECONTAMINATION_WITH_BWA {
    take:
        reads           // tuple(meta, reads)
        ref_genome      // tuple(meta, ref_genome, ref_genome_index)

    main:
    input = ref_genome
    if (!ref_genome.map{it -> it[2]}) {
        INDEX_FASTA(ref_genome)
        input = INDEX_FASTA.out.fasta_with_index
    }
    ALIGNMENT(reads, input)

    SAMTOOLS_FASTQ(ALIGNMENT.out.bams.map{it -> [it[0], it[2]]}, false)

    emit:
        decontaminated_reads = SAMTOOLS_FASTQ.out.fastq
}

process GREP_READS {

    container 'quay.io/biocontainers/pigz:2.3.4'

    tag "${meta.id}"

    input:
    tuple val(meta), path(input_ch_reads)
    tuple val(meta), path(remove_list)

    output:
    tuple val(meta), path("filtered*") , emit: cleaned_reads

    script:

    // define reads
    reads = input_ch_reads.collect()
    def input_reads = "";
    if ( meta.single_end ) {
        """
        zgrep -A 3 -F -v -f ${remove_list} ${reads[0]} | pigz > filtered_${reads[0]}
        """
    } else {
        """
        zgrep -A 3 -F -v -f ${remove_list} ${reads[0]} | pigz > filtered_${reads[0]}
        zgrep -A 3 -F -v -f ${remove_list} ${reads[1]} | pigz > filtered_${reads[1]}
        """
    }

}

workflow DECONTAMINATION_WITH_BMTAGGER {
    take:
        reads           // tuple(meta, reads)
        ref_genome_bitmask
        ref_genome_srprism

    main:
    input_format = channel.value("fastq")
    output_directory = channel.value("bmtagger_output")
    BMTAGGER(reads, ref_genome_bitmask, ref_genome_srprism, input_format, output_directory)
    GREP_READS(reads, BMTAGGER.out.output_host)

    emit:
        decontaminated_reads = GREP_READS.out.cleaned_reads
}
