
process ASSESSMCPPROPORTIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.5--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(fwd_flag), val(rev_flag), path(fastq)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def var_region = "${meta.var_region}"
    def strands = ""

    if (fwd_flag == "auto" && rev_flag == "auto") {
        strands = "FR"
    } else if (fwd_flag == "auto"){
        strands = "F"
    } else if (rev_flag == "auto") {
        strands = "R"
    }

    if (strands == ""){
        """
        touch ${prefix}_${var_region}_mcp_cons.tsv
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mgnify-pipelines-toolkit: 0.1.5
        END_VERSIONS
        """
    } else{
        """
        assess_mcp_proportions \\
            -i ${fastq} \\
            -s ${prefix}_${var_region} \\
            -st ${strands} \\
            -o ./

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mgnify-pipelines-toolkit: 0.1.5
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def var_region = "${meta.var_region}"

    """
    touch ${prefix}_${var_region}_mcp_cons.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: 0.1.5
    END_VERSIONS
    """
}
