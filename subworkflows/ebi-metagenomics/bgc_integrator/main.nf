include { PYRODIGAL                                      } from '../../../modules/ebi-metagenomics/pyrodigal/main'
include { ANTISMASH                                      } from '../../../modules/ebi-metagenomics/antismash/main'
include { ANTISMASH_JSON_TO_GFF                          } from '../../../modules/ebi-metagenomics/antismash_json_to_gff/main'
include { CONCATENATE_GFFS as CONCATENATE_ANTISMASH_GFFS } from '../../../modules/ebi-metagenomics/concatenate_gffs/main'
include { INTERPROSCAN                                   } from '../../../modules/ebi-metagenomics/interproscan/main'
include { SANNTIS                                        } from '../../../modules/ebi-metagenomics/sanntis/main'
include { CONCATENATE_GFFS as CONCATENATE_SANNTIS_GFFS   } from '../../../modules/ebi-metagenomics/concatenate_gffs/main'
include { GECCO_RUN                                      } from '../../../modules/nf-core/gecco/run/main'
include { GECCO_CONVERT                                  } from '../../../modules/nf-core/gecco/convert/main'
include { CONCATENATE_GFFS as CONCATENATE_GECCO_GFFS     } from '../../../modules/ebi-metagenomics/concatenate_gffs/main'
include { MAPPER                                         } from '../../../modules/ebi-metagenomics/bgc_mapper/main'

workflow BGC_INTEGRATOR {

    take:
    ch_inputs           // channel: tuple( val(meta), path(assembly), val(integ_type) )

    main:
    ch_versions = Channel.empty()

    ch_assem = ch_inputs.map{ meta, contigs, _integ_type -> tuple(meta, contigs) }

    PYRODIGAL( ch_assem, 'gbk' )
    ch_versions = ch_versions.mix(PYRODIGAL.out.versions)

    INTERPROSCAN(
        PYRODIGAL.out.faa,
        [file(params.interproscan_database, checkIfExists: true), params.interproscan_database_version],
    )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    SANNTIS(
        INTERPROSCAN.out.ips_annotations.join( PYRODIGAL.out.faa ) 
    )
    ch_versions = ch_versions.mix(SANNTIS.out.versions)

    GECCO( PYRODIGAL.out.annotations, file(params.gecco_hmm, checkIfExists: true) )
    ch_versions = ch_versions.mix(GECCO.out.versions)

    def antismash_channel = channel.empty()
    antismash_channel = ch_assem.join( PYRODIGAL.out.fna ).join( PYRODIGAL.out.annotations )
    ANTISMASH( antismash_channel , antismash_db )
    ch_versions = ch_versions.mix(ANTISMASH.out.versions)

    
    ANTISMASH_JSON_TO_GFF(
        ANTISMASH.out.json
    )
    ch_versions = ch_versions.mix(ANTISMASH_JSON_TO_GFF.out.versions.first())


    MAPPER(
        SANNTIS.out.gff
        .join(GECCO.out.gff)
        .join(ANTISMASH.out.gff)
        .join(integ_type)
    )
    ch_versions = ch_versions.mix(INTEGRATOR.out.versions)


    emit:
    gff      = INTEGRATOR.out.gff           // channel: [ val(meta), [ bam ] ]
    versions = ch_versions                  // channel: [ versions.yml ]
}
