include { DIAMOND_BLASTP  } from '../../../modules/ebi-metagenomics/diamond/blastp/main'
include { CATPACK_CONTIGS } from '../../../modules/ebi-metagenomics/catpack/contigs/main'

workflow CONTIGS_TAXONOMIC_CLASSIFICATION {
    take:
    contigs     // [ val(meta), path(assembly_fasta) ]
    proteins    // [ val(meta), path(proteins_fasta) ]
    cat_db      // [ val(meta), path(catdb_folder)  ]
    taxonomy_db // [ val(meta), path(cattax_folder) ]

    main:

    ch_versions = Channel.empty()

    /*
    * Note:
    * The CAT tool does not use the provided cat_db, as the alignment is performed by the Diamond step.
    * However, this option is mandatory in the CAT source code.
    */

    DIAMOND_BLASTP(
        proteins,
        [[id: "cat-db"], file("${cat_db[1]}/*.dmnd", checkIfExists: true)],
        6, // blast - txt
        []
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions.first())

    CATPACK_CONTIGS(
        contigs,
        cat_db,
        taxonomy_db,
        proteins,
        DIAMOND_BLASTP.out.txt
    )
    ch_versions = ch_versions.mix(CATPACK_CONTIGS.out.versions.first())

    emit:
    diamond_blast_tsv         = DIAMOND_BLASTP.out.txt
    contig2classification_tsv = CATPACK_CONTIGS.out.contig2classification
    versions                  = ch_versions
}
