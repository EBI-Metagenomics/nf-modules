// Subworkflow to generate BGCs annotation using Sanntis, Gecco and Antismash
// Outputs are integrated into a single GFF3 format output
include { GFF2GBK                               } from '../../../modules/ebi-metagenomics/gff2gbk/main'
include { ANTISMASH_ANTISMASH                   } from '../../../modules/nf-core/antismash/antismash/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES  } from '../../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { ANTISMASHGFFBUILDER                   } from '../../../modules/ebi-metagenomics/antismashgffbuilder/main'
include { SANNTIS                               } from '../../../modules/ebi-metagenomics/sanntis/main'
include { GECCO_RUN                             } from '../../../modules/nf-core/gecco/run/main'
include { GECCO_CONVERT                         } from '../../../modules/nf-core/gecco/convert/main'
include { BGCSMAPPER                            } from '../../../modules/ebi-metagenomics/bgcsmapper/main'
include { INTERPROSCAN                          } from '../../../modules/ebi-metagenomics/interproscan/main'
include { ARIA2                                 } from '../../../modules/nf-core/aria2/main'
include { UNTAR                                 } from '../../../modules/nf-core/untar/main'

workflow BGC_INTEGRATOR {

    take:
    ch_inputs           // channel: tuple( val(meta), path(contigs), path(gff), path(proteins), path(ips_annot) )
    ch_antismash_db     // channel: path( antishmash_db )
    ch_ips_db           // channel: path( interproscan_db )
    interproscan_db_url // channel: val( ips_databse_url )
    skip_sanntis        // boolean
    skip_gecco          // boolean
    skip_antismash      // boolean

    main:
    ch_versions = Channel.empty()

    // Extract individual inputs from input channel
    ch_gff = ch_inputs.map{ meta, _contigs, gff, _proteins, _ips_annot -> tuple(meta, gff) }
    ch_togbk_input = ch_inputs.map{ meta, contigs, gff, proteins, _ips_annot -> tuple(meta, contigs, gff, proteins) }
    ch_ips = ch_inputs.map{ meta, _contigs, _gff, _proteins, ips_annot -> tuple(meta, ips_annot) }
    ch_prots = ch_inputs.map{ meta, _contigs, _gff, proteins, _ips_annot -> tuple(meta, proteins) }

    // Transforming gff into gbk format
    GFF2GBK( ch_togbk_input )

    // Running BGCs prediction tools
    // Sanntis
    /*
     * ─────────────────────────────────────────────────────────────────────
     * IPS logic for SanntiS:
     * - If ips_annot is not a path -> run InterProScan to generate it
     * - If ips_annot is provided -> use it
     * ─────────────────────────────────────────────────────────────────────
     */

    if (!skip_sanntis) {

        ch_ips.view { meta, ips -> "IPS debug: ${meta.id} -> ${ips} :: ${ips?.getClass()?.name}" }

        // Split IPS input into "provided" vs "missing" based on whether it's a Path
        ch_ips.branch { meta, ips_annot ->
            provided: ips_annot instanceof Path
            missing:  !(ips_annot instanceof Path)
        }.set { ch_ips_branched }

        // InterProScan input needs proteins only; we join to keep sample keys aligned with "missing"
        ch_prots_missing = ch_prots
            .join(ch_ips_branched.missing)
            .map { meta, prots, _missing -> tuple(meta, prots) }

        // Only download/run InterProScan if there are samples missing IPS annotations
        if (ch_ips_db) {
            INTERPROSCAN( ch_prots_missing, ch_ips_db )
            ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)
            ch_ips_built = INTERPROSCAN.out.tsv
        } else {
            ARIA2( [ [ id:'interproscan_db' ], interproscan_db_url ] )
            ch_versions = ch_versions.mix(ARIA2.out.versions.first())

            UNTAR( ARIA2.out.downloaded_file )
            ch_interproscan_db = UNTAR.out.untar.map { it[1] }

            INTERPROSCAN( ch_prots_missing, ch_interproscan_db )
            ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)
            ch_ips_built = INTERPROSCAN.out.tsv
        }

        // Unified IPS channel for SanntiS: provided + built
        ips_tsv = ch_ips_branched.provided.mix(ch_ips_built)

        // Run SanntiS
        ch_sanntis_input = GFF2GBK.out.gbk
            .join(ips_tsv)
            .map { meta, gbk, ips ->
                [ meta, gbk, ips, [] ]
            }

        SANNTIS( ch_sanntis_input )
        ch_versions = ch_versions.mix(SANNTIS.out.versions)
        ch_sanntis_results = SANNTIS.out.gff
    }

    // Gecco
    if (!skip_gecco) {
        GECCO_RUN( GFF2GBK.out.gbk )
        ch_versions = ch_versions.mix(GECCO_RUN.out.versions)

        GECCO_CONVERT( GECCO_RUN.out )

        ch_gecco_input = GECCO_RUN.out.clusters
            .mix(GECCO_RUN.out.gbk)
            .groupTuple(by:0)
            .map { meta, paths ->
                [meta, paths[0], paths[1]]
            }
        GECCO_CONVERT( ch_gecco_input, "clusters", "gff" )
        ch_gecco_results = GECCO_CONVERT.out.gff
    }

    // AntiSMASH
    if (!skip_antismash && ch_antismash_db) {
        antismash_db = ch_antismash_db
    } else if (!skip_antismash && !ch_antismash_db) {
        ANTISMASH_ANTISMASHDOWNLOADDATABASES()
        ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.versions)
        antismash_db = ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.database
    }

    if (!skip_antismash) {
        ANTISMASH_ANTISMASH( GFF2GBK.out.gbk, antismash_db, [] )
        ch_versions = ch_versions.mix(ANTISMASH_ANTISMASH.out.versions)

        ANTISMASHGFFBUILDER( ANTISMASH_ANTISMASH.out.json )
        ch_antismash_results = ANTISMASHGFFBUILDER.out.gff
    }

    // Gathering samples with results
    // Ensure each optional result channel exists even if tool is skipped
    ch_sanntis_results    = skip_sanntis   ? Channel.empty() : ch_sanntis_results
    ch_gecco_results      = skip_gecco     ? Channel.empty() : ch_gecco_results
    ch_antismash_results  = skip_antismash ? Channel.empty() : ch_antismash_results

    // Per-tool flags: tuple(meta, true)
    ch_sanntis_flags = skip_sanntis
        ? Channel.empty()
        : ch_sanntis_results.map { meta, _file -> tuple(meta, true) }

    ch_gecco_flags = skip_gecco
        ? Channel.empty()
        : ch_gecco_results.map { meta, _file -> tuple(meta, true) }

    ch_antismash_flags = skip_antismash
        ? Channel.empty()
        : ch_antismash_results.map { meta, _file -> tuple(meta, true) }

    // Combine flags -> per-sample has_results boolean
    ch_has_bgc_results = ch_sanntis_flags
        .mix(ch_gecco_flags)
        .mix(ch_antismash_flags)
        .groupTuple()
        .map { meta, _vals -> tuple(meta, true) }

    // Filter base GFF to samples with >=1 tool output
    ch_gff_filtered = ch_gff
        .join(ch_has_bgc_results, remainder: true)
        .filter { meta, gff, has_results -> has_results == true }
        .map    { meta, gff, _has_results -> tuple(meta, gff) }

    /*
    * ─────────────────────────────────────────────────────────────────────
    * Build BGCSMAPPER input with remainder joins (missing tool output -> [])
    * Expected BGCSMAPPER input: tuple(meta, base_gff, sanntis_gff, gecco_gff, antismash_gff)
    * ─────────────────────────────────────────────────────────────────────
    */

    // Start with tuple(meta, base_gff)
    ch_for_bgcsmapper = ch_gff_filtered

    // Add SanntiS (or [])
    if (!skip_sanntis) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_sanntis_results, remainder: true)
            .map { meta, gff, sanntis_gff -> tuple(meta, gff, sanntis_gff ?: []) }
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff -> tuple(meta, gff, []) }
    }

    // Add GECCO (or [])
    if (!skip_gecco) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_gecco_results, remainder: true)
            .map { meta, gff, sanntis_gff, gecco_gff -> tuple(meta, gff, sanntis_gff, gecco_gff ?: []) }
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff, sanntis_gff -> tuple(meta, gff, sanntis_gff, []) }
    }

    // Add antiSMASH (or [])
    if (!skip_antismash) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_antismash_results, remainder: true)
            .map { meta, gff, sanntis_gff, gecco_gff, antismash_gff -> tuple(meta, gff, sanntis_gff, gecco_gff, antismash_gff ?: []) }
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff, sanntis_gff, gecco_gff -> tuple(meta, gff, sanntis_gff, gecco_gff, []) }
    }

    // Run integrator only on samples with >=1 prediction (ch_for_bgcsmapper is already filtered)
    BGCSMAPPER(ch_for_bgcsmapper)

    // Handle cases where no samples pass the filter (integrator produces no output)
    ch_gff_output = BGCSMAPPER.out.gff.ifEmpty([])

    emit:
    gff      = ch_gff_output           // channel: [ val(meta), [ gff ] ]
    versions = ch_versions             // channel: [ versions.yml ]
}
