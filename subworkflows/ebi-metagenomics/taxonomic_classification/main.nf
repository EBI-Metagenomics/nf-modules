include { DIAMOND_BLASTP  } from '../../../modules/ebi-metagenomics/diamond/blastp/main'
include { CATPACK_CONTIGS } from '../../../modules/ebi-metagenomics/catpack/contigs/main'

workflow TAXONOMIC_CLASSIFICATION {
    take:
    contigs                     // [ val(meta), path(file) ]
    proteins                    // [ val(meta), path(file) ]
    cat_db                      // [ val(meta), path(db_folder)  ]
    taxonomy_db                 // [ val(meta), path(tax_folder) ]
 
    main:

    ch_versions = Channel.empty()

    Channel.of(cat_db)
        .map { meta, folder ->
            [meta, file("$folder/*.dmnd")]
        }
        .set {diamond_db}

    DIAMOND_BLASTP(proteins, diamond_db.first(), 6, [])
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions.first())

    // CAT does not use provided cat_db because alignment is already done by Diamond
    // However this option is obligatory in the CAT source python code
    CATPACK_CONTIGS(contigs, cat_db, taxonomy_db, proteins, DIAMOND_BLASTP.out.txt)
    ch_versions = ch_versions.mix(CATPACK_CONTIGS.out.versions.first())

    emit:
    diamond_tsv = DIAMOND_BLASTP.out.txt
    contig2classification_tsv = CATPACK_CONTIGS.out.contig2classification
    versions = ch_versions

}
