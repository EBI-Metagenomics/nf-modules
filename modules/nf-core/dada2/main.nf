process DADA2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0':
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*asvs.fasta"), emit: dada2_out
    path("*asv_counts.tsv")             , emit: dada2_asv_counts
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_args = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    $moduleDir/dada2.R \\
        $args \\
        $prefix \\
        $reads_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -1 | cut -d' ' -f3 )
        dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.asvs.fasta
    touch ${prefix}.asv_counts.tsv

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
    """
}
