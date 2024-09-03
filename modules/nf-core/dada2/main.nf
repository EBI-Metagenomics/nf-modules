// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

process DADA2 {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0':
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*asvs.fasta")    , emit: dada2_out
    path("*asv_counts.tsv")                 , emit: dada2_asv_counts
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( meta.single_end ){
        """
        $moduleDir/dada2.R \\
            $args \\
            $prefix \\
            $reads
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    } else {
        """
        $moduleDir/dada2.R \\
            $args \\
            $prefix \\
            ${reads[0]} \\
            ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    }

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
