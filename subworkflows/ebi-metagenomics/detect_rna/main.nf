// General subwf to detect RNA by using provided RFAM models. It uses cmscan or cmsearch.
// Output: deoverlapped table and chosen fasta file with RNA sequences.

// Use cmscan mode if input fasta file is small and models file is quite big (usecase: mags-catalogues-pipeline)
// Important note: .cm file should be cmpress-ed before execution
// Use cmsearch mode if input fasta is massive and models file contains chosen set of models (usecase: ASA)

/* NF-CORE */
include { SEQKIT_SPLIT2                             } from '../../../modules/nf-core/seqkit/split2/main'
include { CAT_CAT as CONCATENATE_CMSEARCH_DEOVERLAP } from '../../../modules/nf-core/cat/cat/main'

/* EBI-METAGENOMICS */
include { INFERNAL_CMSEARCH                         } from '../../../modules/ebi-metagenomics/infernal/cmsearch/main'
include { INFERNAL_CMSCAN                           } from '../../../modules/ebi-metagenomics/infernal/cmscan/main'
include { CONVERTCMSCANTOCMSEARCH                   } from '../../../modules/ebi-metagenomics/convertcmscantocmsearch/main'
include { CMSEARCHTBLOUTDEOVERLAP                   } from '../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main'
include { EASEL_ESLSFETCH                           } from '../../../modules/ebi-metagenomics/easel/eslsfetch/main'
include { EXTRACTCOORDS                             } from '../../../modules/ebi-metagenomics/extractcoords/main'


workflow DETECT_RNA {

    take:
    ch_fasta          // channel: [ val(meta), [ fasta ] ]
    rfam              // folder: rfam for cmsearch/cmscan
    claninfo          // file: claninfo for cmsearchtbloutdeoverlap
    mode              // cmsearch/cmscan
    separate_subunits // val: boolean
    // analysis_type     // val: str (either AMPLICON or ASSEMBLY)

    main:

    ch_versions = Channel.empty()
    cmsearch_ch = Channel.empty()

    // Chunk the fasta into files with at most params.proteins_chunksize sequences
    SEQKIT_SPLIT2(
        ch_fasta
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    ch_sequences = SEQKIT_SPLIT2.out.reads.transpose()

    if ( mode == 'cmsearch' ) {
        INFERNAL_CMSEARCH(
            ch_sequences,
            rfam
        )
        ch_versions = ch_versions.mix(INFERNAL_CMSEARCH.out.versions.first())
        cmsearch_ch = INFERNAL_CMSEARCH.out.cmsearch_tbl
    }
    else if (mode == 'cmscan') {
       INFERNAL_CMSCAN(
            ch_sequences,
            rfam
       )
       ch_versions = ch_versions.mix(INFERNAL_CMSCAN.out.versions.first())

       CONVERTCMSCANTOCMSEARCH(INFERNAL_CMSCAN.out.cmscan_tbl)
       ch_versions = ch_versions.mix(CONVERTCMSCANTOCMSEARCH.out.versions.first())

       cmsearch_ch = CONVERTCMSCANTOCMSEARCH.out.cmsearch_tblout
    }

    CMSEARCHTBLOUTDEOVERLAP(
        cmsearch_ch,
        claninfo
    )
    ch_versions = ch_versions.mix(CMSEARCHTBLOUTDEOVERLAP.out.versions.first())

    CONCATENATE_CMSEARCH_DEOVERLAP(
        CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_CMSEARCH_DEOVERLAP.out.versions.first())

    ch_easel = ch_fasta
                .join(CONCATENATE_CMSEARCH_DEOVERLAP.out.file_out)
    EASEL_ESLSFETCH(
        ch_easel
    )
    ch_versions = ch_versions.mix(EASEL_ESLSFETCH.out.versions.first())

    EXTRACTCOORDS(
        EASEL_ESLSFETCH.out.easel_coords,
        EASEL_ESLSFETCH.out.matched_seqs_with_coords,
        separate_subunits
    )
    ch_versions = ch_versions.mix(EXTRACTCOORDS.out.versions.first())

    emit:
    cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped   // channel: [ val(meta), [ deoverlapped ] ]
    easel_out              = EASEL_ESLSFETCH.out.easel_coords                           // channel: [ val(meta), [ fasta ] ]
    ssu_fasta              = EXTRACTCOORDS.out.ssu_fasta                                // channel: [ val(meta), [ fasta ] ]
    lsu_fasta              = EXTRACTCOORDS.out.lsu_fasta                                // channel: [ val(meta), [ fasta ] ]
    concat_ssu_lsu_coords  = EXTRACTCOORDS.out.concat_ssu_lsu_coords                    // channel: [ val(meta), [ txt ] ]
    versions               = ch_versions                                                // channel: [ versions.yml ]
}

