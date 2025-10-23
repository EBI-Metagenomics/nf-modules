 /*
  Subworkflow to annotate eukaryotic genes
*/
include {REPEATMODELER_BUILDDATABASE} from '../../../modules/ebi-metagenomics/repeatmodeler/builddatabase/main.nf'
include {REPEATMODELER_REPEATMODELER} from '../../../modules/ebi-metagenomics/repeatmodeler/repeatmodeler/main.nf'
include {REPEATMASKER_REPEATMASKER} from '../../../modules/ebi-metagenomics/repeatmasker/repeatmasker/main.nf'
include {BRAKER3} from '../../../modules/ebi-metagenomics/braker3/main.nf'
include {METAEUK_EASYPREDICT} from '../../../modules/ebi-metagenomics/metaeuk/easypredict/main.nf'
include {EUKGENEOVERLAP} from '../../../modules/ebi-metagenomics/eukgeneoverlap/main.nf'

workflow EUK_GENE_CALLING {
    take:
        genome_info // tuple(meta, genome: fasta)
        protein_evidence
        transcriptomic_evidence
        overlap_threshold
        

    main:
        // de novo repeat identification and masking
        REPEATMODELER_BUILDDATABASE(
            genome_info
        )

        REPEATMODELER_REPEATMODELER(
            REPEATMODELER_BUILDDATABASE.out.db
        )

        REPEATMASKER_REPEATMASKER(
            genome_info,
            REPEATMODELER_REPEATMODELER.out.fasta
        )


        // call genes based on input evidence
        if ( protein_evidence && transcriptomic_evidence ) {
            BRAKER3(
                REPEATMASKER_REPEATMASKER.out.masked,
                null,
                transcriptomic_evidence,
                null,
                protein_evidence,
                null
            )
            METAEUK_EASYPREDICT(
                REPEATMASKER_REPEATMASKER.out.masked,
                null
            )
            EUKGENEOVERLAP(
                genome_info.meta,
                BRAKER3.out.gff3,
                METAEUK_EASYPREDICT.out.gff,
                BRAKER3.out.aa,
                METAEUK_EASYPREDICT.out.faa,
                overlap_threshold
            )
        }

        if ( protein_evidence && !transcriptomic_evidence ) {
            BRAKER3(
                REPEATMASKER_REPEATMASKER.out.masked,
                null,
                null,
                null,
                protein_evidence,
                null
            )

            METAEUK_EASYPREDICT(
                REPEATMASKER_REPEATMASKER.out.masked,
                null
            )
            EUKGENEOVERLAP(
                genome_info.meta,
                BRAKER3.out.gff3,
                METAEUK_EASYPREDICT.out.gff,
                BRAKER3.out.aa,
                METAEUK_EASYPREDICT.out.faa,
                overlap_threshold
            )
            }

        if ( !protein_evidence && transcriptomic_evidence ) {
            BRAKER3(
                REPEATMASKER_REPEATMASKER.out.masked,
                null,
                transcriptomic_evidence,
                null,
                null,
                null
            )
        }

        if ( !protein_evidence && !transcriptomic_evidence ) {
            BRAKER3(
                REPEATMASKER_REPEATMASKER.out.masked,
                null,
                null,
                null,
                null,
                null
            )
        }


    emit:
        braker_gff = BRAKER3.out.gff3
        braker_proteins = BRAKER3.out.aa
        braker_ffn = BRAKER3.out.cds
        metaeuk_gff = METAEUK_EASYPREDICT.out.gff
        metaeuk_proteins = METAEUK_EASYPREDICT.out.faa
        metaeuk_ffn = METAEUK_EASYPREDICT.out.codon
        softmasked_genome = REPEATMASKER_REPEATMASKER.out.masked
}