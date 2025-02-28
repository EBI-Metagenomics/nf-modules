process BRAKER_BRAKER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/braker3:3.0.8--hdfd78af_0':
        'biocontainers/braker3:3.0.8--hdfd78af_0' }"

    input:
    tuple val(meta), path(genome)
    path bam
    path fa
    path gff

    output:
    tuple val(meta), path("${prefix}/*.gtf"), emit: gtf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def accepted_hits = bam ? "--bam ${bam}" : ''
    def proteins = fa ? "--prot_seq ${fa}" : ''
    def hints = gff ? "--hints ${gff}" : ''
    """
    braker.pl \\
        $args \\
        --genome ${genome} \\
        $accepted_hits \\
        $proteins \\
        $hints \\
        --threads $task.cpus \\
        --workingdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker: \$(echo \$(braker.pl --version) | sed -E 's/[^0-9]*([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    touch $prefix/${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker: \$(echo \$(braker.pl --version) | sed -E 's/[^0-9]*([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
