
include { SAMTOOLS_FAIDX                       } from '../../../modules/ebi-metagenomics/samtools/faidx/main'
include { TABIX_BGZIP as TABIX_BGZIP_FA        } from '../../../modules/ebi-metagenomics/tabix/bgzip/main'

workflow INDEX_FASTA {

    take:
    ch_fasta  // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    TABIX_BGZIP_FA(ch_fasta)
    ch_versions = ch_versions.mix(TABIX_BGZIP_FA.out.versions.first())

    SAMTOOLS_FAIDX(TABIX_BGZIP_FA.out.output)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    emit:
    fasta_gz       = TABIX_BGZIP_FA.out.output      // channel: [ val(meta), [ gz ] ]
    fa_gz_gzi      = TABIX_BGZIP_FA.out.gzi         // channel: [ val(meta), [ fai ] ]
    fai            = SAMTOOLS_FAIDX.out.fai         // channel: [ val(meta), [ fai ] ]
    versions = ch_versions                          // channel: [ versions.yml ]
}
