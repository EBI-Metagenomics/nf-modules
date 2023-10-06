process DEEPTMHMM {
    label 'process_medium'

    conda "bioconda::pybiolib=1.1.1342"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybiolib:1.1.1342--pyhdfd78af_0':
        'biocontainers/pybiolib:1.1.1342--pyhdfd78af_0' }"

    input:
    path fasta

    output:
    path "biolib_results/deeptmhmm_results.md"      , emit: md
    path "biolib_results/predicted_topologies.3line", emit: line3
    path "biolib_results/TMRs.gff3"                 , emit: gff3
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""

    """
    biolib \\
        run \\
        DTU/DeepTMHMM \\
        --fasta ${fasta} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biolib: \$(echo \$(biolib --version) | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""

    """
    mkdir biolib_results
    touch biolib_results/deeptmhmm_results.md
    touch biolib_results/predicted_topologies.3line
    touch biolib_results/TMRs.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biolib: \$(echo \$(biolib --version) | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p' )
    END_VERSIONS
    """
}
