// Subworkflow to generate antimicrobial resistance annotation from protein and gene files
// Outputs are standardised and integrated into a single GFF3 format output

/* NF-CORE */

include { AMRFINDERPLUS_UPDATE             } from '../../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                } from '../../../modules/nf-core/amrfinderplus/run/main'
include { DEEPARG_DOWNLOADDATA             } from '../../../modules/nf-core/deeparg/downloaddata/main'
include { DEEPARG_PREDICT                  } from '../../../modules/nf-core/deeparg/predict/main'
include { RGI_CARDANNOTATION               } from '../../../modules/nf-core/rgi/cardannotation/main'
include { RGI_MAIN                         } from '../../../modules/nf-core/rgi/main/main'
include { UNTAR as UNTAR_CARD              } from '../../../modules/nf-core/untar/main'
include { HAMRONIZATION_RGI                } from '../../../modules/nf-core/hamronization/rgi/main'
include { HAMRONIZATION_DEEPARG            } from '../../../modules/nf-core/hamronization/deeparg/main'
include { AMR_INTEGRATOR                   } from '../../../modules/nf-core/amrintegrator/main'

workflow AMR_ANNOTATION {
    take:
    ch_inputs            // channel: tuple( val(meta), path(aminoacids), path(cds_gff) )

    main:
    ch_versions = Channel.empty()

    // Extract individual components from input channel
    ch_faa = ch_inputs.map{ meta, aminoacids, _cds_gff -> tuple(meta, aminoacids) }
    ch_gff = ch_inputs.map{ meta, _aminoacids, cds_gff -> tuple(meta, cds_gff) }


    // RUNNING ANNOTATION TOOLS
    // AMRfinderplus run
    // Prepare channel for database
    if (!params.arg_skip_amrfinderplus && params.arg_amrfinderplus_db) {
        ch_amrfinderplus_db = Channel
            .fromPath(params.arg_amrfinderplus_db, checkIfExists: true)
            .first()
    }
    else if (!params.arg_skip_amrfinderplus && !params.arg_amrfinderplus_db) {
        AMRFINDERPLUS_UPDATE()
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    if (!params.arg_skip_amrfinderplus) {
        AMRFINDERPLUS_RUN(ch_faa, ch_amrfinderplus_db)
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)
    }


    // RGI run
    if (!params.arg_skip_rgi) {
        if (!params.arg_rgi_db) {

            // Download and untar CARD
            UNTAR_CARD([[], file('https://card.mcmaster.ca/latest/data', checkIfExists: true)])
            ch_versions = ch_versions.mix(UNTAR_CARD.out.versions)
            rgi_db = UNTAR_CARD.out.untar.map { it[1] }
            RGI_CARDANNOTATION(rgi_db)
            card = RGI_CARDANNOTATION.out.db
            ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)
        }
        else {

            // Use user-supplied database
            rgi_db = file(params.arg_rgi_db, checkIfExists: true)
            if (!rgi_db.contains("card_database_processed")) {
                RGI_CARDANNOTATION(rgi_db)
                card = RGI_CARDANNOTATION.out.db
                ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)
            }
            else {
                card = rgi_db
            }
        }

        RGI_MAIN(ch_faa, card, [])
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

        // Reporting
        HAMRONIZATION_RGI(RGI_MAIN.out.tsv, 'tsv', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions)
    }

    // DeepARG prepare download
    if (!params.arg_skip_deeparg && params.arg_deeparg_db) {
        ch_deeparg_db = Channel
            .fromPath(params.arg_deeparg_db, checkIfExists: true)
            .first()
    }
    else if (!params.arg_skip_deeparg && !params.arg_deeparg_db) {
        DEEPARG_DOWNLOADDATA()
        ch_versions = ch_versions.mix(DEEPARG_DOWNLOADDATA.out.versions)
        ch_deeparg_db = DEEPARG_DOWNLOADDATA.out.db
    }

    // DeepARG run
    if (!params.arg_skip_deeparg) {

        ch_faa
            .map { it ->
                def meta = it[0]
                def anno = it[1]
                def model = params.arg_deeparg_model

                [meta, anno, model]
            }
            .set { ch_input_for_deeparg }

        DEEPARG_PREDICT(ch_input_for_deeparg, ch_deeparg_db)
        ch_versions = ch_versions.mix(DEEPARG_PREDICT.out.versions)

        // Reporting
        // Note: currently hardcoding versions as unreported by DeepARG
        // Make sure to update on version bump.
        ch_input_to_hamronization_deeparg = DEEPARG_PREDICT.out.arg.mix(DEEPARG_PREDICT.out.potential_arg)
        HAMRONIZATION_DEEPARG(ch_input_to_hamronization_deeparg, 'tsv', '1.0.4', params.arg_deeparg_db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions)
    }

    // Integrate and transform into a single GFF3 output file
    AMR_INTEGRATOR(
            HAMRONIZATION_DEEPARG.out.tsv
        .join(
            HAMRONIZATION_RGI.out.tsv
        ).join(
            AMRFINDERPLUS_RUN.out.tsv
        ).join(
            ch_gff
    )

    emit:
    gff      = AMRINTEGRATOR.out.gff           // channel: [ val(meta), [ gff ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
