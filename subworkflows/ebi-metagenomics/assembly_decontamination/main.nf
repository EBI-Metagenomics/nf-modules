include { BLAST_BLASTN } from '../../../modules/ebi-metagenomics/blast/blastn/main'
include { SEQKIT_GREP } from '../../../modules/ebi-metagenomics/seqkit/grep/main'

workflow ASSEMBLY_DECONTAMINATION {

    take:
    assembly               // [ val(meta), path(assembly_fasta) ]
    reference_genome       // [ path(reference_genome) ]

    main:
    ch_versions = Channel.empty()
    cleaned_contigs = Channel.empty()
    ch_blast_ref = Channel.fromPath( "${params.blast_reference_genomes_folder}/${reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": reference_genome], files ]
            }

    BLAST_BLASTN(assembly,ch_blast_ref)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    SEQKIT_GREP(assembly,BLAST_BLASTN.out.txt.flatten().last())
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    cleaned_contigs = SEQKIT_GREP.out.filter
    versions        = ch_versions 
}
