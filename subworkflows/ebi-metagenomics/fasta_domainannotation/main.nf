include { BLAST_MAKEBLASTDB } from '../../../modules/ebi-metagenomics/blast/makeblastdb/main'
include { BLAST_BLASTP      } from '../../../modules/ebi-metagenomics/blast/blastp/main'
include { INTERPROSCAN      } from '../../../modules/ebi-metagenomics/interproscan/main'
include { EGGNOGMAPPER      } from '../../../modules/ebi-metagenomics/eggnogmapper/main'

workflow FASTA_DOMAINANNOTATION {

    take:
    ch_fasta       // channel: [ val(meta), [ fasta ] ]
    ch_blast_fasta // channel: [ blast_fasta ]
    ch_eggnog      // channel: [ eggnog_db, eggnog_data_dir, diamond_db ]

    main:

    ch_versions = Channel.empty()

    BLAST_MAKEBLASTDB ( blast_fasta )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions.first())
    BLAST_BLASTP ( ch_fasta, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_BLASTP.out.versions.first())

    out_ext = 'tsv'
    INTERPROSCAN ( ch_fasta, out_ext )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.first())

    EGGNOGMAPPER ( ch_fasta, ch_eggnog.eggnog_db, ch_eggnog.eggnog_data_dir, ch_eggnog.diamond_db )
    ch_versions = ch_versions.mix(EGGNOGMAPPER.out.versions.first())

    emit:
    blastp_csv       = BLAST_BLASTP.out.csv // channel: [ val(meta), [ csv ] ]
    inteproscan_tsv  = INTERPROSCAN.out.tsv // channel: [ val(meta), [ tsv ] ]
    eggnogmapper_csv = EGGNOGMAPPER.out.csv // channel: [ val(meta), [ csv ] ]

    versions = ch_versions // channel: [ versions.yml ]
}

