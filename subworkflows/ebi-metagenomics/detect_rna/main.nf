
include { INFERNAL_CMSEARCH           } from '../../../modules/ebi-metagenomics/infernal/cmsearch/main'
include { INFERNAL_CMSCAN             } from '../../../modules/ebi-metagenomics/infernal/cmscan/main'
include { CONVERSTCMSCANTOCMSEARCH    } from '../../../modules/ebi-metagenomics/convertcmscantocmsearch/main'
include { CMSEARCHTBLOUTDEOVERLAP     } from '../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main'
include { EASEL_ESLSFETCH             } from '../../../modules/ebi-metagenomics/easel/eslsfetch/main'

workflow DETECT_RNA {

    take:
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    rfam         // folder: rfam for cmsearch/cmscan
    claninfo     // file: claninfo for cmsearchtbloutdeoverlap
    mode         // cmsearch/cmscan

    main:

    ch_versions = Channel.empty()
    cmsearch_ch = Channel.empty()

    if ( mode == 'cmsearch' ) {
        INFERNAL_CMSEARCH(
            ch_fasta,
            rfam
        )
        ch_versions = ch_versions.mix(INFERNAL_CMSEARCH.out.versions.first())
        cmsearch_ch = INFERNAL_CMSEARCH.out.cmsearch_tbl
    }
    else if (mode == 'cmscan') {
       INFERNAL_CMSCAN(
            ch_fasta,
            rfam
       )
       ch_versions = ch_versions.mix(INFERNAL_CMSCAN.out.versions.first())

       CONVERSTCMSCANTOCMSEARCH(INFERNAL_CMSCAN.out.cmscan_tbl)
       // todo add version of toolkit
       cmsearch_ch = CONVERSTCMSCANTOCMSEARCH.out.cmsearch_tblout
    }

    CMSEARCHTBLOUTDEOVERLAP(
        cmsearch_ch,
        claninfo
    )
    ch_versions = ch_versions.mix(CMSEARCHTBLOUTDEOVERLAP.out.versions.first())

    ch_easel = ch_fasta
                .join(CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped)
    EASEL_ESLSFETCH(
        ch_easel
    )
    ch_versions = ch_versions.mix(EASEL_ESLSFETCH.out.versions.first())

    emit:
    cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped   // channel: [ val(meta), [ deoverlapped ] ]
    easel_out              = EASEL_ESLSFETCH.out.easel_coords                           // channel: [ val(meta), [ fasta ] ]
    versions               = ch_versions                                                // channel: [ versions.yml ]
}

