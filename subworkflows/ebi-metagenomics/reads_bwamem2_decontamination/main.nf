include { BWAMEM2_MEM     } from '../../../modules/ebi-metagenomics/bwamem2/mem/main'
include { BWAMEM2_INDEX   } from '../../../modules/ebi-metagenomics/bwamem2/index/main'
include { SAMTOOLS_BAM2FQ } from '../../../modules/ebi-metagenomics/samtools/bam2fq/main'


workflow READS_BWAMEM2_DECONTAMINATION   {

    take:
    ch_reads       // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_reference   // channel (mandatory): [ val(meta2), path(ref_index) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()

    // Checking if the bwamem2 index files are present
    def expected_index_files = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]

    def index_files_count = expected_index_files.count { pattern ->
        ch_reference.map { meta2, ref_index -> ref_index.count { file -> file.endsWith(pattern) } }.sum()
    }    

    println index_files_count
    ch_reference.view()

    if (index_files_count == expected_index_files.size()) {
        BWAMEM2_MEM(ch_reads, ch_reference)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    } else {
        reference_fasta = ch_reference.map { meta2, ref_index -> [ meta2, ref_index[0]] }
        BWAMEM2_INDEX(reference_fasta)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        BWAMEM2_MEM(ch_reads, BWAMEM2_INDEX.out.index)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    }

    SAMTOOLS_BAM2FQ(BWAMEM2_MEM.out.bam.map { meta, bam, bai -> [ meta, bam ] }, ch_reads.map { meta, reads -> meta.single_end == false } )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())
   
    emit:
    decontaminated_reads = SAMTOOLS_BAM2FQ.out.reads  // channel: [ val(meta), [ path(decont_reads) ]]
    versions = ch_versions                            // channel: [ versions.yml ]

}



