process CREATE_DB_BMTAGGER {

    label 'process_low'

    conda "bioconda::bmtagger==3.101"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bmtagger:3.101--h470a237_4':
        'biocontainers/bmtagger:3.101--h470a237_4' }"

    input:
    path(reference_fasta)

    output:
    path("${reference_fasta.baseName}.bitmask")         , emit: bitmask
    path("${reference_fasta.baseName}.srprism.*")       , emit: srprism
    path("${reference_fasta}.n*")              , emit: blast_db

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bmtool -d ${reference_fasta} -o ${reference_fasta.baseName}.bitmask -A 0 -w 18
    srprism mkindex -i ${reference_fasta} -o ${reference_fasta.baseName}.srprism -M 7168
    makeblastdb -in ${reference_fasta} -dbtype nucl
    """

    stub:
    """
    touch reference.bitmask reference.srprism.idx reference.fasta.nsq
    """
}
