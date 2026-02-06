// Subworkflow to generate BGCs annotation using Sanntis, Gecco and Antismash
// Outputs are integrated into a single GFF3 format output
include { GFF2GBK                                        } from '../../../modules/ebi-metagenomics/gff2gbk/main'
include { ANTISMASH_ANTISMASH                            } from '../../../modules/nf-core/antismash/antismash/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES           } from '../../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { ANTISMASH_JSON_TO_GFF                          } from '../../../modules/ebi-metagenomics/antismash_json_to_gff/main'
include { SANNTIS                                        } from '../../../modules/ebi-metagenomics/sanntis/main'
include { GECCO_RUN                                      } from '../../../modules/nf-core/gecco/run/main'
include { GECCO_CONVERT                                  } from '../../../modules/nf-core/gecco/convert/main'
include { BGCSINTEGRATOR                                 } from '../../../modules/ebi-metagenomics/bgcsintegrator/main'

workflow BGC_INTEGRATOR {

    take:
    ch_inputs           // channel: tuple( val(meta), path(contigs), path(gff), path(proteins), path(ips_annot) )
    ch_antishmash_db    // channel: path( antishmash_db )

    main:
    ch_versions = Channel.empty()

    // Extract individual inputs from input channel
    ch_togbk_input = ch_inputs.map{ meta, contigs, gff, proteins, _ips_annot -> tuple(meta, contigs, gff, proteins) }
    ch_ips = ch_inputs.map{ meta, _contigs, _gff, proteins, ips_annot -> tuple(meta, ips_annot) }

    // Preparing databases
    if (ch_antishmash_db) {
        antishmash_db = ch_antishmash_db
    } else {
        ANTISMASH_ANTISMASHDOWNLOADDATABASES()
        ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.versions)
        antishmash_db = ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.database
    }

    // Transforming gff into gbk format
    GFF2GBK( ch_togbk_input )

    // Running BGCs prediction tools
    ch_sanntis_input = GFF2GBK.out.gbk
        .join(ch_ips)
        .map { meta, gbk, ips ->
            [meta, gbk, ips, []]
        }
    SANNTIS( ch_sanntis_input )
    ch_versions = ch_versions.mix(SANNTIS.out.versions)

    GECCO_RUN( GFF2GBK.out.gbk )
    ch_versions = ch_versions.mix(GECCO_RUN.out.versions)

    ch_gecco_conv_input = GECCO_RUN.out.clusters
        .mix(GECCO_RUN.out.gbk)
        .groupTuple(by:0)
        .map { meta, paths ->
            [meta, paths[0], paths[1]]
        }
    GECCO_CONVERT( ch_gecco_conv_input, "clusters", "gff" )
    ch_versions = ch_versions.mix(GECCO_CONVERT.out.versions)

    ANTISMASH_ANTISMASH( GFF2GBK.out.gbk, antismash_db, [] )
    ch_versions = ch_versions.mix(ANTISMASH_ANTISMASH.out.versions)

    ANTISMASH_JSON_TO_GFF( ANTISMASH_ANTISMASH.out.json )
    ch_versions = ch_versions.mix(ANTISMASH_JSON_TO_GFF.out.versions.first())

    BGCSINTEGRATOR(
        SANNTIS.out.gff
        .join(GECCO_CONVERT.out.gff)
        .join(ANTISMASH_JSON_TO_GFF.out.gff)
    )
    ch_versions = ch_versions.mix(BGCSINTEGRATOR.out.versions)

    // Handle cases where no predictions are made (integrator produces no output)
    ch_gff_output = BGCSINTEGRATOR.out.gff.ifEmpty([])

    emit:
    gff      = ch_gff_output           // channel: [ val(meta), [ bam ] ]
    versions = ch_versions             // channel: [ versions.yml ]
}
