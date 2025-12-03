// Subworkflow to generate antimicrobial resistance annotation from protein and gene files
// Outputs are standardised and integrated into a single GFF3 format output
include { AMRFINDERPLUS_UPDATE             } from '../../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                } from '../../../modules/nf-core/amrfinderplus/run/main'
include { WGET as WGET_DEEPARG             } from '../../../modules/nf-core/wget/main'
include { UNZIP as UNZIP_DEEPARG           } from '../../../modules/nf-core/unzip/main'
include { DEEPARG_PREDICT                  } from '../../../modules/nf-core/deeparg/predict/main'
include { RGI_CARDANNOTATION               } from '../../../modules/nf-core/rgi/cardannotation/main'
include { RGI_MAIN                         } from '../../../modules/nf-core/rgi/main/main'
include { WGET as WGET_CARD                } from '../../../modules/nf-core/wget/main'
include { UNTAR as UNTAR_CARD              } from '../../../modules/nf-core/untar/main'
include { HAMRONIZATION_RGI                } from '../../../modules/nf-core/hamronization/rgi/main'
include { HAMRONIZATION_DEEPARG            } from '../../../modules/nf-core/hamronization/deeparg/main'
include { AMRINTEGRATOR                    } from '../../../modules/ebi-metagenomics/amrintegrator/main'

workflow AMR_ANNOTATION {
    take:
    ch_inputs                  // channel: tuple( val(meta), path(aminoacids), path(cds_gff) )
    arg_amrfinderplus_db       // channel: path( amrfinderplus_db )
    arg_deeparg_db             // channel: path( deeparg_db )
    arg_deeparg_db_version     // channel: val( deeparg_db_version )
    arg_deeparg_model          // channel: val( deeparg_model )
    arg_rgi_db                 // channel: path( rgi_db )
    arg_skip_amrfinderplus     // boolean 
    arg_skip_deeparg           // boolean
    arg_skip_rgi               // boolean

    main:
    ch_versions = channel.empty()

    // Extract individual components from input channel
    ch_faa = ch_inputs.map{ meta, aminoacids, _cds_gff ->
        // Add is_proteins flag to meta
        def meta_with_protein_flag = meta + [is_proteins: true]
        tuple(meta_with_protein_flag, aminoacids)
    }
    ch_gff = ch_inputs.map{ meta, _aminoacids, cds_gff -> tuple(meta, cds_gff) }

    // RUNNING ANNOTATION TOOLS
    // AMRfinderplus run
    // Prepare channel for database
    if (!arg_skip_amrfinderplus && arg_amrfinderplus_db) {
        ch_amrfinderplus_db = arg_amrfinderplus_db
    }
    else if (!arg_skip_amrfinderplus && !arg_amrfinderplus_db) {
        AMRFINDERPLUS_UPDATE()
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions.first())
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    if (!arg_skip_amrfinderplus) {
        AMRFINDERPLUS_RUN(ch_faa, ch_amrfinderplus_db)
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions.first())
        ch_amrfinderplus_results = AMRFINDERPLUS_RUN.out.report
    }

    // RGI run
    if (!arg_skip_rgi) {
        if (!arg_rgi_db) {
            // Download and untar CARD
            ch_card_url = channel.of([
                [id: 'card_database'],
                'https://card.mcmaster.ca/latest/data'
            ])
            WGET_CARD(ch_card_url)
            ch_versions = ch_versions.mix(WGET_CARD.out.versions.first())
            UNTAR_CARD(WGET_CARD.out.outfile)
            ch_versions = ch_versions.mix(UNTAR_CARD.out.versions.first())
            rgi_db = UNTAR_CARD.out.untar.map { it[1] }
            RGI_CARDANNOTATION(rgi_db)
            card = RGI_CARDANNOTATION.out.db
            ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions.first())
        }
        else {
            // Use user-supplied database
            rgi_db = arg_rgi_db
            if (!rgi_db.contains("card_database_processed")) {
                RGI_CARDANNOTATION(rgi_db)
                card = RGI_CARDANNOTATION.out.db
                ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions.first())
            }
            else {
                card = rgi_db
            }
        }

        RGI_MAIN(ch_faa, card, [])
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions.first())

        // Reporting
        HAMRONIZATION_RGI(RGI_MAIN.out.tsv, 'tsv', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions.first())
        ch_rgi_results = HAMRONIZATION_RGI.out.tsv
    }

    // DeepARG prepare download
    if (!arg_skip_deeparg && arg_deeparg_db) {
        ch_deeparg_db = arg_deeparg_db
    }
    else if (!arg_skip_deeparg && !arg_deeparg_db) {
        ch_deeparg_url = channel.of([
            [ id: 'deeparg_db' ],
            'https://zenodo.org/records/8280582/files/deeparg.zip?download=1'
        ])
    
        // Download the database using WGET
        WGET_DEEPARG( ch_deeparg_url )
        ch_versions = ch_versions.mix(WGET_DEEPARG.out.versions.first())
    
        // Unzip the downloaded database
        UNZIP_DEEPARG(WGET_DEEPARG.out.outfile )
        ch_versions = ch_versions.mix(UNZIP_DEEPARG.out.versions.first())
    
        // Extract the unzipped directory path
        ch_deeparg_db = UNZIP_DEEPARG.out.unzipped_archive.map { meta, dir -> dir }
    }

    // DeepARG run
    if (!arg_skip_deeparg) {
        ch_faa
            .map { it ->
                def meta = it[0]
                def anno = it[1]
                def model = arg_deeparg_model
                [meta, anno, model]
            }
            .set { ch_input_for_deeparg }

        DEEPARG_PREDICT(ch_input_for_deeparg, ch_deeparg_db)
        ch_versions = ch_versions.mix(DEEPARG_PREDICT.out.versions.first())

        // Reporting
        // Note: currently hardcoding versions as unreported by DeepARG
        // Make sure to update on version bump
        HAMRONIZATION_DEEPARG(DEEPARG_PREDICT.out.arg, 'tsv', '1.0.4', arg_deeparg_db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions.first())
        ch_deeparg_results = HAMRONIZATION_DEEPARG.out.tsv
    }

    // Integrate and transform into a single GFF3 output file
    // Key insight: Tool results have meta with is_proteins:true, but ch_gff has original meta
    // We need to join by ID to avoid meta mismatch

    // Create keyed channels using meta.id for consistent joining
    ch_gff_keyed = ch_gff.map { meta, gff -> [meta.id, meta, gff] }

    // Create flag channels for tools that ran (keyed by meta.id)
    ch_deeparg_flags = arg_skip_deeparg ? 
        channel.empty() : 
        ch_deeparg_results.map { meta, _file -> [meta.id, 'deeparg'] }
    
    ch_rgi_flags = arg_skip_rgi ? 
        channel.empty() : 
        ch_rgi_results.map { meta, _file -> [meta.id, 'rgi'] }
    
    ch_amrfinder_flags = arg_skip_amrfinderplus ? 
        channel.empty() : 
        ch_amrfinderplus_results.map { meta, _file -> [meta.id, 'amrfinder'] }

    // Combine all flags to identify samples that have at least one tool result
    ch_all_flags = ch_deeparg_flags
        .mix(ch_rgi_flags)
        .mix(ch_amrfinder_flags)
        .groupTuple()
        .map { id, tools -> [id, true] }  // Create a simple boolean flag

    // Filter GFF channel to only samples with results
    ch_gff_filtered = ch_gff_keyed
        .join(ch_all_flags, remainder: true)
        .filter { _id, _meta, _gff, has_results ->
            has_results == true
        }
        .map { id, meta, gff, _has_results ->
            [id, meta, gff]
        }

    // Build the input channel for AMRINTEGRATOR
    // Join tool results (keyed by meta.id) with the filtered GFF
    
    // Add deeparg results (or empty list if skipped)
    if (!arg_skip_deeparg) {
        ch_deeparg_keyed = ch_deeparg_results.map { meta, file -> [meta.id, file] }
        ch_for_amrintegrator = ch_gff_filtered
            .join(ch_deeparg_keyed, remainder: true)
            .map { id, meta, gff, deeparg -> 
                [id, meta, gff, deeparg ?: []]  // Replace null with empty list
            }
    } else {
        ch_for_amrintegrator = ch_gff_filtered
            .map { id, meta, gff -> [id, meta, gff, []] }  // Empty list for skipped tool
    }

    // Add rgi results (or empty list if skipped)
    if (!arg_skip_rgi) {
        ch_rgi_keyed = ch_rgi_results.map { meta, file -> [meta.id, file] }
        ch_for_amrintegrator = ch_for_amrintegrator
            .join(ch_rgi_keyed, remainder: true)
            .map { id, meta, gff, deeparg, rgi -> 
                [id, meta, gff, deeparg, rgi ?: []]  // Replace null with empty list
            }
    } else {
        ch_for_amrintegrator = ch_for_amrintegrator
            .map { id, meta, gff, deeparg -> 
                [id, meta, gff, deeparg, []]  // Empty list for skipped tool
            }
    }

    // Add amrfinderplus results (or empty list if skipped)
    if (!arg_skip_amrfinderplus) {
        ch_amrfinder_keyed = ch_amrfinderplus_results.map { meta, file -> [meta.id, file] }
        ch_for_amrintegrator = ch_for_amrintegrator
            .join(ch_amrfinder_keyed, remainder: true)
            .map { id, meta, gff, deeparg, rgi, amrfinder -> 
                [id, meta, gff, deeparg, rgi, amrfinder ?: []]  // Replace null with empty list
            }
    } else {
        ch_for_amrintegrator = ch_for_amrintegrator
            .map { id, meta, gff, deeparg, rgi -> 
                [id, meta, gff, deeparg, rgi, []]  // Empty list for skipped tool
            }
    }

    // Final reordering to match AMRINTEGRATOR input: [meta, deeparg, rgi, amrfinder, gff]
    // Remove the id key and reorder
    ch_for_amrintegrator = ch_for_amrintegrator
        .map { id, meta, gff, deeparg, rgi, amrfinder ->
            [meta, deeparg, rgi, amrfinder, gff]
        }

    AMRINTEGRATOR(ch_for_amrintegrator)
    ch_versions = ch_versions.mix(AMRINTEGRATOR.out.versions.first())

    emit:
    gff      = AMRINTEGRATOR.out.gff            // channel: [ val(meta), [ gff ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
