include { BLAST_MAKEBLASTDB } from '../../../modules/ebi-metagenomics/blast/makeblastdb/main'
include { BLAST_BLASTP      } from '../../../modules/ebi-metagenomics/blast/blastp/main'
include { DIAMOND_MAKEDB    } from '../../../modules/ebi-metagenomics/diamond/makedb/main'
include { DIAMOND_BLASTP    } from '../../../modules/ebi-metagenomics/diamond/blastp/main'
include { INTERPROSCAN      } from '../../../modules/ebi-metagenomics/interproscan/main'

workflow FASTA_DOMAINANNOTATION {

    take:
    ch_fasta       // channel: [ val(meta), path(fasta) ]
    ch_blast_fasta // channel: /path/to/reference/fasta for blast
    val_blast_mode // value: blast or diamond

    main:

    ch_versions = Channel.empty()

    if (val_blast_mode == "blast") {
        println("blast mode")
        BLAST_MAKEBLASTDB ( ch_blast_fasta )
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
        BLAST_BLASTP ( ch_fasta, BLAST_MAKEBLASTDB.out.db )
        ch_versions = ch_versions.mix(BLAST_BLASTP.out.versions)
        blastp_csv = BLAST_BLASTP.out.csv
    } else if (val_blast_mode == "diamond") {
        println("diamond mode")
        DIAMOND_MAKEDB ( ch_blast_fasta )
        ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
        out_ext = 'txt'
        blast_columns = '' // defaults: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
        DIAMOND_BLASTP ( ch_fasta, DIAMOND_MAKEDB.out.db, out_ext, blast_columns )
        ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)
        blastp_csv = DIAMOND_BLASTP.out.txt
    } else {
        throw new Exception("Invalid mode value '$val_blast_mode'. Should be 'blast' or 'diamond'.")
    }

    out_ext = 'tsv'
    INTERPROSCAN ( ch_fasta, out_ext )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    emit:
    blastp_csv       = blastp_csv // channel: [ val(meta), [ csv ] ]
    inteproscan_tsv  = INTERPROSCAN.out.tsv // channel: [ val(meta), [ tsv ] ]
    versions         = ch_versions // channel: [ versions.yml ]
}

