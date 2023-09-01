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
    def args_bmtool = task.ext.args_bmtool ?: '-A 0 -w 18'
    def args_srprism = task.ext.args_srprism ?: '-M 7168'
    def args_makeblastdb = task.ext.args_makeblastdb ?: '-dbtype nucl'

    script:
    """
    bmtool -d ${reference_fasta} -o ${reference_fasta.baseName}.bitmask ${args_bmtool}
    srprism mkindex -i ${reference_fasta} -o ${reference_fasta.baseName}.srprism ${args_srprism}
    makeblastdb -in ${reference_fasta} ${args_makeblastdb}

    bmtool_version=\$(bmtool --version | tr ' ' '\\t' | cut -f2)
    srprism_version="2.4.24-alpha"  //\$(srprism --version | grep version | tr ' ' '\\t' | cut -f2)
    blast_version=\$(makeblastdb -h | grep -w "Application" | tr ',' '\\t' | cut -f2 | tr ' ' '\\t' | cut -f3)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bmtool: \${bmtool_version}
        srprism: \${srprism_version}
        makeblastdb: \${blast_version}
    END_VERSIONS
    """

    stub:
    """
    touch reference.bitmask reference.srprism.idx reference.fasta.nsq

    bmtool_version=\$(bmtool --version | tr ' ' '\\t' | cut -f2)
    srprism_version="2.4.24-alpha"  //\$(srprism --version | grep version | tr ' ' '\\t' | cut -f2)
    blast_version=\$(makeblastdb -h | grep -w "Application" | tr ',' '\\t' | cut -f2 | tr ' ' '\\t' | cut -f3)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bmtool: \${bmtool_version}
        srprism: \${srprism_version}
        makeblastdb: \${blast_version}
    END_VERSIONS
    """
}
