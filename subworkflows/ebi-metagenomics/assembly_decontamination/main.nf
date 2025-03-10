include { BLAST_BLASTN } from '../../../modules/ebi-metagenomics/blast/blastn/main'
include { SEQKIT_GREP  } from '../../../modules/ebi-metagenomics/seqkit/grep/main'

workflow ASSEMBLY_DECONTAMINATION {

    take:
    assembly           // [ val(meta), path(assembly_fasta) ]
    ch_blast_ref       // [ val(meta2), path(reference_genome_files) ]

    main:
    ch_versions = Channel.empty()

    BLAST_BLASTN(assembly,ch_blast_ref)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    SEQKIT_GREP(assembly,BLAST_BLASTN.out.txt.flatten().last())
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    cleaned_contigs = SEQKIT_GREP.out.filter
    versions        = ch_versions 
}
