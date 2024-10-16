include { SAMTOOLS_FAIDX } from '../../../../modules/ebi-metagenomics/samtools/faidx/main'
include { TABIX_BGZIP    } from '../../../../modules/ebi-metagenomics/tabix/bgzip/main'

workflow INDEX_FASTA {

    take:
    ch_fasta  // channel: [ val(meta), [ fasta ] ]
    output_dir   // The directory path to save output files

    main:

    ch_versions = Channel.empty()

    // Run bgzip on fasta files
    TABIX_BGZIP(ch_fasta)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    // Index the gzipped fasta files
    SAMTOOLS_FAIDX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    // Extract paths from the tuple before passing them to the publish process
    fasta_gz       = TABIX_BGZIP.out.output.map { it[1] }.flatten()        // channel: [ val(meta), [ gz ] ]
    fai            = SAMTOOLS_FAIDX.out.fai.map { it[1] }.flatten()        // channel: [ val(meta), [ fai ] ]
    fa_gz_gzi      = SAMTOOLS_FAIDX.out.gzi.map { it[1] }.flatten()        // channel: [ val(meta), [ gzi ] ]
    versions       = ch_versions                                           // channel: [ versions.yml ]

    // Call the process to publish files
    PUBLISH_OUTPUT_FILES(fasta_gz, fai, fa_gz_gzi, output_dir)
}

// PUBLISH_OUTPUT_FILES process to save the output files
process PUBLISH_OUTPUT_FILES {

    input:
    path fasta_gz
    path fa_gz_gzi
    path fai
    val output_dir

    script:
    """
    echo "Files in working directory before copying:"
    ls -lh

    mkdir -p ${output_dir}

    echo "Copying files..."
    cp $fasta_gz ${output_dir}/
    cp $fa_gz_gzi ${output_dir}/
    cp $fai ${output_dir}/

    """
}

