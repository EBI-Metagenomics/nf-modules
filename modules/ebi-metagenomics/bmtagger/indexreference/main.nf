process BMTAGGER_INDEXREFERENCE {
    tag "${reference_fasta}"
    label 'process_low'

    conda "bioconda::bmtagger=3.101"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bmtagger:3.101--h470a237_4':
        'biocontainers/bmtagger:3.101--h470a237_4' }"

    input:
    path(reference_fasta)

    output:
    path("${reference_fasta.simpleName}.bitmask")  , emit: bitmask
    path("${reference_fasta.simpleName}.srprism.*"), emit: srprism
    path("${reference_fasta.simpleName}.n*")       , emit: blast_db
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args_bmtool = task.ext.args_bmtool ?: '-A 0 -w 18'
    def args_srprism = task.ext.args_srprism ?: '-M 7168'
    def args_makeblastdb = task.ext.args_makeblastdb ?: '-dbtype nucl'
    """
    bmtool -d ${reference_fasta} -o ${reference_fasta.simpleName}.bitmask ${args_bmtool}
    srprism mkindex -i ${reference_fasta} -o ${reference_fasta.simpleName}.srprism ${args_srprism}
    makeblastdb -in ${reference_fasta} ${args_makeblastdb}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bmtool: \$(bmtool --version | tr ' ' '\\t' | cut -f2)
        srprism: \$(srprism 2>&1 | grep version | sed 's/^.*version //; s/ .*\$//')
        makeblastdb: \$(makeblastdb -version 2>&1 | sed 's/^.*makeblastdb: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${reference_fasta.simpleName}.bitmask
    touch ${reference_fasta.simpleName}.srprism.idx
    touch ${reference_fasta.simpleName}.srprism.amp
    touch ${reference_fasta.simpleName}.srprism.idx
    touch ${reference_fasta.simpleName}.srprism.imp
    touch ${reference_fasta.simpleName}.srprism.map
    touch ${reference_fasta.simpleName}.srprism.pmp
    touch ${reference_fasta.simpleName}.srprism.rmp
    touch ${reference_fasta.simpleName}.srprism.ss
    touch ${reference_fasta.simpleName}.srprism.ssa
    touch ${reference_fasta.simpleName}.srprism.ssd
    touch ${reference_fasta.simpleName}.nsq
    touch ${reference_fasta.simpleName}.nhr
    touch ${reference_fasta.simpleName}.nhr
    touch ${reference_fasta.simpleName}.nin
    touch ${reference_fasta.simpleName}.nin
    touch ${reference_fasta.simpleName}.nsq
    touch ${reference_fasta.simpleName}.nsq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bmtool: \$(bmtool --version | tr ' ' '\\t' | cut -f2)
        srprism: \$(srprism 2>&1 | grep version | sed 's/^.*version //; s/ .*\$//')
        makeblastdb: \$(makeblastdb -version 2>&1 | sed 's/^.*makeblastdb: //; s/ .*\$//')
    END_VERSIONS
    """
}
