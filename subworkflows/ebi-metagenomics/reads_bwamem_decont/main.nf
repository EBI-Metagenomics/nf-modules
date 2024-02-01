// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BWAMEM2_MEM        } from '../../../modules/ebi-metagenomics/bwamem2/mem'
include { SAMTOOLS_FASTQ     } from '../../../modules/nf-core/samtools/fastq/main'


workflow READS_BWAMEM_DECONT {

    take:
        reads               // tuple(meta, reads)
        ref_genome
        ref_genome_index

    main:
    ch_versions = Channel.empty()

    to_align = reads.map { meta, reads -> 
        [ meta, reads, ref_genome, ref_genome_index ]
    }

    BWAMEM2_MEM( to_align )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    SAMTOOLS_BAM2FQ( ALIGNMENT.out.bam.map { meta, ref_fasta, bam, bai -> [ meta, bam ] }, true )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    emit:
    decontaminated_reads = SAMTOOLS_BAM2FQ.out.fastq  // channel: [ val(meta), [ decont_reads ]]
    versions = ch_versions                            // channel: [ versions.yml ]

}

