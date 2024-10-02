
include { SEQFU_CHECK            } from '../../../modules/ebi-metagenomics/seqfu/check/main'
include { ASSESSMCPPROPORTIONS   } from '../../../modules/ebi-metagenomics/assessmcpproportions/main'
include { FASTQSUFFIXHEADERCHECK } from '../../../modules/ebi-metagenomics/fastqsuffixheadercheck/main'
include { FASTP                  } from '../../../modules/ebi-metagenomics/fastp/main'
include { SEQTK_SEQ              } from '../../../modules/ebi-metagenomics/seqtk/seq/main'

workflow  READS_QC {

    take:
    ch_reads    // channel: [ val(meta), [ fastq ] ]
    save_merged // channel:  val(boolean)

    main:
    ch_versions = Channel.empty()
    
    // ***** Necessary mapping functions *****
    filterBySeqFuStatus = { meta, seqfu_res ->
            seqfu_check_status = seqfu_res[0]
            if (seqfu_check_status == "OK"){
                [ meta ]
            }            
         }

    filterBySuffixHeaderStatus = { meta, sufhd_res ->
            if (sufhd_res.countLines() == 0){
                [ meta ]
            }
        }

    add_mcp_flags = { meta, fastq ->
            [ meta, "auto", "auto", fastq ]
        }

    prepareForMCPCheck = { meta, fwd_flag, rev_flag, fastq ->
            if (meta.single_end){
                [ meta, fwd_flag, rev_flag, fastq ]
            }
            else{
                 [ meta, fwd_flag, rev_flag, fastq[0] ]
            }
        }

    filterByAmpliconCheck = { meta, strategy ->
            if (strategy == "AMPLICON"){
                [ meta ]
            }
        }
    // ***** Necessary mapping functions *****

    SEQFU_CHECK(ch_reads)

    passed_seqfu_reads = SEQFU_CHECK.out.tsv
                   .splitCsv(sep: "\t", elem: 1)
                   .map(filterBySeqFuStatus)
                   .join(ch_reads)

    FASTQSUFFIXHEADERCHECK(passed_seqfu_reads)

    passed_suffixheader_reads = FASTQSUFFIXHEADERCHECK.out.json
                    .map(filterBySuffixHeaderStatus)
                    .join(ch_reads)

    assess_mcp_proportions_input = passed_suffixheader_reads
                    .map(add_mcp_flags)
                    .map(prepareForMCPCheck)

    ASSESSMCPPROPORTIONS(assess_mcp_proportions_input, true)

    passed_amplicon_check_reads = ASSESSMCPPROPORTIONS.out.library_check_out
                    .map(filterByAmpliconCheck)
                    .join(ch_reads)

    FASTP ( passed_amplicon_check_reads, params.save_trimmed_fail, save_merged )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_se_fastp_reads = FASTP
                        .out.reads
                        .filter { it[0].single_end }

    ch_reads_se_and_merged = ch_se_fastp_reads
                            .mix(FASTP.out.reads_merged)

    SEQTK_SEQ(ch_reads_se_and_merged)
    ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions.first())

    emit:
    reads               = FASTP.out.reads           // channel: [ val(meta), [ fastq ] ]
    reads_se_and_merged = ch_reads_se_and_merged    // channel: [ val(meta), [ fastq ] ]
    fastp_summary_json  = FASTP.out.json            // channel: [ val(meta), [ json ] ]
    reads_fasta         = SEQTK_SEQ.out.fastx       // channel: [ val(meta), [ fasta ] ]
    versions            = ch_versions               // channel: [ versions.yml ]
}

