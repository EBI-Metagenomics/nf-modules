/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Decontamination subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM2_MEM        } from '../../../modules/ebi-metagenomics/bwamem2/mem'
include { SAMTOOLS_BAM2FQ     } from '../../../modules/ebi-metagenomics/samtools/bam2fq/main'


workflow READS_BWAMEM_DECONT {

    take:
    ch_reads       // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index       // channel (mandatory): [ path(index) ]

    main:
    ch_versions = Channel.empty()

    BWAMEM2_MEM( ch_reads, ch_index )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    def split_arg = "";
    if ( ch_reads.collect().size() == 2 ) {
        split_arg = TRUE;
    } else {
        split_arg = FALSE;
    }

    SAMTOOLS_BAM2FQ( BWAMEM2_MEM.out.bam, BWAMEM2_MEM.out.bai, split_arg )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    emit:
    decontaminated_reads = SAMTOOLS_BAM2FQ.out.reads  // channel: [ val(meta), [ decont_reads ]]
    versions = ch_versions                            // channel: [ versions.yml ]

}

