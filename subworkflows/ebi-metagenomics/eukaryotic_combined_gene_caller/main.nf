 /*
  Subworkflow to annotate eukaryotic genes
*/
include {REPEATMODELER_BUILDDATABASE} from '../../../modules/ebi-metagenomics/repeatmodeler/builddatabase/main.nf'
include {REPEATMODELER_REPEATMODELER} from '../../../modules/ebi-metagenomics/repeatmodeler/repeatmodeler/main.nf'
include {REPEATMASKER_REPEATMASKER} from '../../../modules/ebi-metagenomics/repeatmasker/repeatmasker/main.nf'
include {BRAKER3} from '../../../modules/ebi-metagenomics/braker3/main.nf'
include {METAEUK_EASYPREDICT} from '../../../modules/ebi-metagenomics/metaeuk/easypredict/main.nf'
include {EUKGENEOVERLAP} from '../../../modules/ebi-metagenomics/eukgeneoverlap/main.nf'
include {EUKGENERENAME} from '../../../modules/ebi-metagenomics/eukgenerename/main.nf'

workflow EUK_GENE_CALLING {
    take:
        genome_info // tuple(meta, genome: fasta)
        protein_evidence
        transcriptomic_evidence
        overlap_threshold
        

    main:

        // mask repeats in the genome
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


        overlap_results = Channel.empty()
        braker_overlap_proteins = Channel.empty()
        metaeuk_overlap_proteins = Channel.empty()
        braker_unique_proteins = Channel.empty()
        metaeuk_unique_proteins = Channel.empty()

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
            
            braker_gff = BRAKER3.out.gff3
            metaeuk_gff = METAEUK_EASYPREDICT.out.gff
            overlap_results = EUKGENEOVERLAP.out.overlap_tsv
            braker_overlap_proteins = EUKGENEOVERLAP.out.braker_overlap_faa
            braker_overlap_nucelotide = EUKGENEOVERLAP.out.braker_overlap_ffn
            metaeuk_overlap_proteins = EUKGENEOVERLAP.out.metaeuk_overlap_faa
            metaeuk_overlap_nucleotide = EUKGENEOVERLAP.out.metaeuk_overlap_ffn
            braker_proteins = EUKGENEOVERLAP.out.braker_unique_faa
            braker_ffn = EUKGENEOVERLAP.out.braker_unique_ffn
            metaeuk_proteins = EUKGENEOVERLAP.out.metaeuk_unique_faa
            metaeuk_ffn = EUKGENEOVERLAP.out.metaeuk_unique_ffn 
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
            
            braker_gff = BRAKER3.out.gff3
            metaeuk_gff = METAEUK_EASYPREDICT.out.gff
            overlap_results = EUKGENEOVERLAP.out.overlap_tsv
            braker_overlap_proteins = EUKGENEOVERLAP.out.braker_overlap_faa
            braker_overlap_nucelotide = EUKGENEOVERLAP.out.braker_overlap_ffn
            metaeuk_overlap_proteins = EUKGENEOVERLAP.out.metaeuk_overlap_faa
            metaeuk_overlap_nucleotide = EUKGENEOVERLAP.out.metaeuk_overlap_ffn
            braker_proteins = EUKGENEOVERLAP.out.braker_unique_faa
            braker_ffn = EUKGENEOVERLAP.out.braker_unique_ffn
            metaeuk_proteins = EUKGENEOVERLAP.out.metaeuk_unique_faa
            metaeuk_ffn = EUKGENEOVERLAP.out.metaeuk_unique_ffn 
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

            EUKGENERENAME(
                genome_info.meta,
                BRAKER3.out.aa,
                BRAKER3.out.cds
            )

            braker_gff = BRAKER3.out.gff3
            braker_proteins = EUKGENERENAME.out.renamed_braker_aa
            braker_ffn = EUKGENERENAME.out.renamed_braker_ffn
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

            EUKGENERENAME(
                genome_info.meta,
                BRAKER3.out.aa,
                BRAKER3.out.cds
            )

            braker_gff = BRAKER3.out.gff3
            braker_proteins = EUKGENERENAME.out.renamed_braker_aa
            braker_ffn = EUKGENERENAME.out.renamed_braker_ffn
        }


    emit:
        braker_gff = braker_gff
        braker_proteins = braker_proteins
        braker_ffn = braker_ffn
        metaeuk_gff = metaeuk_gff
        metaeuk_proteins = metaeuk_proteins
        metaeuk_ffn = metaeuk_ffn
        overlap_results = overlap_results
        braker_overlap_proteins = braker_overlap_proteins
        metaeuk_overlap_proteins = metaeuk_overlap_proteins
        braker_unique_proteins = braker_unique_proteins
        metaeuk_unique_proteins = metaeuk_unique_proteins
        softmasked_genome = REPEATMASKER_REPEATMASKER.out.masked
}