// Subworkflow to generate antimicrobial resistance annotation from protein and gene files
// Outputs are standardised and integrated into a single GFF3 format output

/* NF-CORE */
include { DEEPARG                     } from '../../../modules/nf-core/deeparg/predict/main'
include { HAMRONIZATION_DEEPARG       } from '../../../modules/nf-core/hamronization/deeparg/main'
include { RGI                         } from '../../../modules/nf-core/rgi/main'
include { HAMRONIZATION_RGI           } from '../../../modules/nf-core/hamronization/rgi/main'
include { AMRFINDERPLUS               } from '../../../modules/nf-core/amrfinderplus/run/main'
include { HAMRONIZATION_AMRFINDERPLUS } from '../../../modules/nf-core/hamronization/amrfinderplus/main'

/* EBI-METAGENOMICS */
include { AMRINTEGRATOR               } from '../../../modules/ebi-metagenomics/amrintegrator/main'


workflow AMR_ANNOTATION {
    ch_inputs            // channel: tuple( val(meta), path(genes), path(aminoacids), path(cds_gff) )
    ch_deeparg_inputs    // channel: tuple( path(deeparg_model), deeparg_version, deeparg_db_version )
    ch_rgi_inputs        // channel: tuple( path(card_db), path(wildcard_db) )
    ch_amrfinderplus_db  // channel: path( amrfinderplus_db)

    main:
    ch_versions = Channel.empty()

    // Extract individual components from input channels
    ch_fnn = ch_inputs.map{ meta, genes, _aminoacids, _cds_gff -> tuple(meta, genes) }
    ch_faa = ch_inputs.map{ meta, _genes, aminoacids, _cds_gff -> tuple(meta, aminoacids) }
    ch_gff = ch_inputs.map{ meta, _genes, _aminoacids, cds_gff -> tuple(meta, cds_gff) }

    // Extract deeparg components
    ch_deeparg_model = ch_deeparg_inputs.map{ deeparg_model, _version, _db_version -> deeparg_model }
    ch_deeparg_version = ch_deeparg_inputs.map{ _deeparg_model, version, _db_version -> version }
    ch_deeparg_db_version = ch_deeparg_inputs.map{ _deeparg_model, _version, db_version -> db_version }

    // Extract rgi components
    ch_rgi_card = ch_rgi_inputs.amp{ card_db, _wildcard_db -> card_db }
    ch_rgi_wildcard = ch_rgi_inputs.amp{ _card_db, wildcard_db -> wildcard_db }


    // RUNNING ANNOTATION TOOLS
    // Run DEEPARG and HAMRONIZATION
    DEEPARG( 
        ch_fnn, 
        'LS', 
        ch_deeparg_model
    )
    ch_versions = ch_versions.mix(DEEPARG.out.versions.first())

    HAMRONIZATION_DEEPARG(
        DEEPARG.out.arg, 
        'tsv', 
        ch_deeparg_version, 
        ch_deeparg_db_version,
    )
    ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions.first())


    // Run RGI and HAMRONIZATION
    RGI( 
        ch_faa,
        ch_rgi_card,
        ch_rgi_wildcard,
    )
    ch_versions = ch_versions.mix(RGI.out.versions.first())

    HAMRONIZATION_RGI(
        RGI.out.tsv,
        'tsv',
        RGI.out.tool_version,
        RGI.out.db_version,
    )
    ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions.first())


    // Run AMRFINDERPLUS and HAMRONIZATION 
    AMRFINDERPLUS(
        ch_faa, 
        ch_amrfinderplus_db,
    )
    ch_versions = ch_versions.mix(AMRFINDERPLUS.out.versions.first())

    HAMRONIZATION_AMRFINDERPLUS(
        AMRFINDERPLUS.out.report,
        'tsv',
        AMRFINDERPLUS.out.tool_version 
        AMRFINDERPLUS.out.db_version,   
    )
    ch_versions = ch_versions.mix(HAMRONIZATION_AMRFINDERPLUS.out.versions.first())
    

    // Integrate and transform into a single GFF3 output file
    AMRINTEGRATOR(
            HAMRONIZATION_DEEPARG.out.tsv
        .join(
            HAMRONIZATION_RGI.out.tsv
        ).join(
            HAMRONIZATION_AMRFINDERPLUS.out.tsv
        ).join(
            ch_gff
    )

    emit:
    //gff      = AMRINTEGRATOR.out.gff           // channel: [ val(meta), [ gff ] ]

    dea = HAMRONIZATION_DEEPARG.out.tsv
    rgi = HAMRONIZATION_RGI.out.tsv
    amr = HAMRONIZATION_AMRFINDERPLUS.out.tsv

    versions = ch_versions                     // channel: [ versions.yml ]
}
