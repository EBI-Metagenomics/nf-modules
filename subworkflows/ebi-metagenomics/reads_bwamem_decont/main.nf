/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Decontamination subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM2_MEM        } from '../../../modules/ebi-metagenomics/bwamem2/mem'
include { SAMTOOLS_BAM2FQ     } from '../../../modules/ebi-metagenomics/samtools/bam2fq/main'


workflow READS_BWAMEM_DECONT {

    take:
        reads               // tuple(meta, reads)
        ref_genome          // path(reference_genome)
        ref_genome_index    // path(reference_genome_index

    main:
    ch_versions = Channel.empty()

    to_align = reads.map { meta, reads -> 
        [ meta, reads, ref_genome, ref_genome_index ]
    }

    BWAMEM2_MEM( to_align )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    SAMTOOLS_BAM2FQ( BWAMEM2_MEM.out.bam.map { meta, bam, bai -> [ meta, bam, bai ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    emit:
    decontaminated_reads = SAMTOOLS_BAM2FQ.out.reads  // channel: [ val(meta), [ decont_reads ]]
    versions = ch_versions                           // channel: [ versions.yml ]

}

