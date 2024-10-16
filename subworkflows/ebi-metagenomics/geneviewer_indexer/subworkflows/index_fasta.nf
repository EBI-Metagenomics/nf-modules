
include { SAMTOOLS_FAIDX    } from '../../../../modules/ebi-metagenomics/samtools/faidx/main'
include { TABIX_BGZIP       } from '../../../../modules/ebi-metagenomics/tabix/bgzip/main'

workflow INDEX_FASTA {

    take:
    ch_fasta  // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()
    ch_fasta.view()

    TABIX_BGZIP(ch_fasta)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    SAMTOOLS_FAIDX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    emit:
    fasta_gz       = TABIX_BGZIP.out.output      // channel: [ val(meta), [ gz ] ]
    fa_gz_gzi      = TABIX_BGZIP.out.gzi         // channel: [ val(meta), [ fai ] ]
    fai            = SAMTOOLS_FAIDX.out.fai         // channel: [ val(meta), [ fai ] ]
    versions = ch_versions                          // channel: [ versions.yml ]
}
