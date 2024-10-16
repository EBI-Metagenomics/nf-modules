include { GFF3_TRIM_FASTA      } from '../../../../modules/ebi-metagenomics/gff3trimfasta/main'
include { SORT_GFF             } from '../../../../modules/ebi-metagenomics/jbrowse/sortgff/main'
include { TABIX_TABIX          } from '../../../../modules/ebi-metagenomics/tabix/tabix/main'
include { TABIX_BGZIP    } from '../../../../modules/ebi-metagenomics/tabix/bgzip/main'

workflow INDEX_GFF {

    take:
        ch_gff  // channel: [ val(meta), [ gff ] ]
        output_dir   // The directory path to save output files

    main:

        ch_versions = Channel.empty()

        GFF3_TRIM_FASTA(ch_gff)
        ch_versions = ch_versions.mix(GFF3_TRIM_FASTA.out.versions.first())

        // Check that GFF3_TRIM_FASTA emits the correct file
        GFF3_TRIM_FASTA.out.gff.view { "Trimmed GFF: $it" }

        SORT_GFF(GFF3_TRIM_FASTA.out.gff)
        ch_versions = ch_versions.mix(SORT_GFF.out.versions.first())

        // Check that SORT_GFF emits the correct sorted file
        SORT_GFF.out.gff.view { "Sorted GFF: $it" }

        TABIX_BGZIP(SORT_GFF.out.gff)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

        TABIX_TABIX(TABIX_BGZIP.out.output)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        // Check that TABIX_TABIX emits the correct tbi files
        TABIX_TABIX.out.tbi.view { "Tabix TBI: $it" }

        gff_gz      = TABIX_BGZIP.out.output.map { it[1] }.flatten()        // channel: [ val(meta), [ gz ] ]
        tbi_files   = TABIX_TABIX.out.tbi.map { it[1] }.flatten()           // Channel: [ val(meta), path(tbi_files) ]
        versions    = ch_versions                                           // Channel: [ versions.yml ]

    // Call the process to publish files
    PUBLISH_OUTPUT_FILES(gff_gz, tbi_files, output_dir)
}

// PUBLISH_OUTPUT_FILES process to save the output files
process PUBLISH_OUTPUT_FILES {

    input:
    path gff_gz
    path tbi_files
    val output_dir

    script:
    """
    echo "Files in working directory before copying:"
    ls -lh

    mkdir -p ${output_dir}

    echo "Copying files..."
    cp $gff_gz ${output_dir}/
    cp $tbi_files ${output_dir}/

    """
}
