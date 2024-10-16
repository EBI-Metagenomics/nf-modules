include { SAMTOOLS_FAIDX } from '../../../../modules/ebi-metagenomics/samtools/faidx/main'
include { TABIX_BGZIP    } from '../../../../modules/ebi-metagenomics/tabix/bgzip/main'

workflow INDEX_FASTA {

    take:
    ch_fasta  // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()
    ch_fasta.view { "Input FASTA: $it" }

    // Run bgzip on fasta files
    TABIX_BGZIP(ch_fasta)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())
    TABIX_BGZIP.out.output.view { "TABIX_BGZIP output: $it" }

    // Index the gzipped fasta files
    SAMTOOLS_FAIDX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
    SAMTOOLS_FAIDX.out.fai.view { "SAMTOOLS fai output: $it" }
    SAMTOOLS_FAIDX.out.gzi.view { "SAMTOOLS gzi output: $it" }

    // Extract paths from the tuple before passing them to the publish process
    fasta_gz       = TABIX_BGZIP.out.output.map { it[1] }.flatten()        // channel: [ val(meta), [ gz ] ]
    fai            = SAMTOOLS_FAIDX.out.fai.map { it[1] }.flatten()        // channel: [ val(meta), [ fai ] ]
    fa_gz_gzi      = SAMTOOLS_FAIDX.out.gzi.map { it[1] }.flatten()        // channel: [ val(meta), [ gzi ] ]
    versions       = ch_versions                                           // channel: [ versions.yml ]

    // Call the process to publish files
    PUBLISH_OUTPUT_FILES(fasta_gz, fai, fa_gz_gzi)
}

// PUBLISH_OUTPUT_FILES process to save the output files
process PUBLISH_OUTPUT_FILES {

    input:
    path fasta_gz
    path fa_gz_gzi
    path fai

    script:
    """
    echo "Files in working directory before copying:"
    ls -lh

    mkdir -p /Users/vikasg/Documents/ws/work/nf-modules/subworkflows/ebi-metagenomics/geneviewer_indexer/output

    echo "Copying files..."
    cp $fasta_gz /Users/vikasg/Documents/ws/work/nf-modules/subworkflows/ebi-metagenomics/geneviewer_indexer/output/
    cp $fa_gz_gzi /Users/vikasg/Documents/ws/work/nf-modules/subworkflows/ebi-metagenomics/geneviewer_indexer/output/
    cp $fai /Users/vikasg/Documents/ws/work/nf-modules/subworkflows/ebi-metagenomics/geneviewer_indexer/output/

    echo "Files in output folder after copying:"
    ls -lh /Users/vikasg/Documents/ws/work/nf-modules/subworkflows/ebi-metagenomics/geneviewer_indexer/output/
    """
}

