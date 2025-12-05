// Subworkflow to generate toxins and virulence factors annotation from protein sequences
// Outputs are filtered by threshold and integrated into a single GFF3 format output
include { PATHOFACT2_DOWNLOADDATA } from '../../../modules/ebi-metagenomics/pathofact2/downloaddata/main'
include { PATHOFACT2_TOXINS       } from '../../../modules/ebi-metagenomics/pathofact2/toxins/main'
include { PATHOFACT2_VIRULENCE    } from '../../../modules/ebi-metagenomics/pathofact2/virulence/main'
include { PATHOFACT2_INTEGRATOR   } from '../../../modules/ebi-metagenomics/pathofact2/integrator/main'

workflow PATHOFACT2 {
    take:
    ch_inputs          // channel: tuple( val(meta), path(aminoacids), path(cds_gff) )
    ch_models          // channel: path( pathofact2_db )

    main:
    ch_versions = channel.empty()

    // Extract individual components from input channel
    ch_faa = ch_inputs.map{ meta, aminoacids, _cds_gff -> tuple(meta, aminoacids) }
    ch_gff = ch_inputs.map{ meta, _aminoacids, cds_gff -> tuple(meta, cds_gff) }

    // Preparing databases
    if (ch_models) {
        pathofact_models = ch_models
    }
    else {
        PATHOFACT2_DOWNLOADDATA()
        ch_versions = ch_versions.mix(PATHOFACT2_DOWNLOADDATA.out.versions.first())
        pathofact_models = PATHOFACT2_DOWNLOADDATA.out.db
    }

    // Running annotation
    PATHOFACT2_TOXINS( ch_faa, pathofact_models )
    ch_versions = ch_versions.mix(PATHOFACT2_TOXINS.out.versions.first())
    
    PATHOFACT2_VIRULENCE( ch_faa, pathofact_models )
    ch_versions = ch_versions.mix(PATHOFACT2_VIRULENCE.out.versions.first())

    // Integrating results in a single gff file
    ch_for_integrator =  ch_gff
        .join(PATHOFACT2_TOXINS.out.tsv)
        .join(PATHOFACT2_VIRULENCE.out.tsv)

    PATHOFACT2_INTEGRATOR( ch_for_integrator )
    ch_versions = ch_versions.mix(PATHOFACT2_INTEGRATOR.out.versions.first())

    emit:
    gff  =  PATHOFACT2_INTEGRATOR.out.gff         // channel: tuple( val(meta), path(gff) )
}
