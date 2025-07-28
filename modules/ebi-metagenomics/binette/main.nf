process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/binette:1.1.2--d469f6785354dd73'
        : 'biocontainers/binette:1.1.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(input, stageAs: "input_bins/*")
    tuple val(meta2), path(contigs)
    val input_type
    path proteins
    path db

    output:
    tuple val(meta), path("final_bins/*")                  , emit: refined_bins
    tuple val(meta), path("final_bins_quality_reports.tsv"), emit: refined_bins_report
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def protein_arg = proteins ? "--proteins ${proteins}" : ""
    def input_arg   = ""

    if (input_type == 'fasta') {
        input_arg = "--bin_dirs input_bins/*"
    } else if (input_type == 'tsv') {
        input_arg = "--contig2bin_tables input_bins/*"
    } else {
        error "Invalid input_type: ${input_type}. Must be 'fasta' or 'tsv'"
    }

    """
    binette \\
        ${input_arg} \\
        ${protein_arg} \\
        --contigs ${contigs} \\
        --threads ${task.cpus} \\
        --checkm2_db ${db} \\
        -o . \\
        -v \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir final_bins/
    touch final_bins/bin_1.fasta
    touch final_bins_quality_reports.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """
}