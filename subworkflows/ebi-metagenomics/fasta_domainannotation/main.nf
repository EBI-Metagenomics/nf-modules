include { BLAST_MAKEBLASTDB } from '../../../modules/ebi-metagenomics/blast/makeblastdb/main'
include { BLAST_BLASTP      } from '../../../modules/ebi-metagenomics/blast/blastp/main'
include { INTERPROSCAN      } from '../../../modules/ebi-metagenomics/interproscan/main'

workflow FASTA_DOMAINANNOTATION {

    take:
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    val_blast_fasta     // channel: /path/to/reference/fasta for blast

    main:

    ch_versions = Channel.empty()

    BLAST_MAKEBLASTDB ( val_blast_fasta )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
    BLAST_BLASTP ( ch_fasta, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_BLASTP.out.versions)

    out_ext = 'tsv'
    INTERPROSCAN ( ch_fasta, out_ext )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    emit:
    blastp_csv       = BLAST_BLASTP.out.csv // channel: [ val(meta), [ csv ] ]
    inteproscan_tsv  = INTERPROSCAN.out.tsv // channel: [ val(meta), [ tsv ] ]

    versions = ch_versions // channel: [ versions.yml ]
}

