process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/binette:1.0.4--pyh7e72e81_0'
        : 'biocontainers/binette:1.1.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(input, stageAs: "input/*")
    tuple val(meta2), path(contigs)
    val(input_type)
    path(proteins)
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
    if (input_type == 'fasta') {
        input_arg = "--bin_dirs input/*"
    } else if (input_type == 'tsv') {
        input_arg = "--contig2bin_tables input/*"
    } else {
        error "Either FASTA directories or TSV files must be provided"
    }
    """
    binette \\
        ${input_arg} \\
        ${protein_arg} \\
        --contigs ${contigs} \\
        -o . \\
        -v \\
        --threads ${task.cpus} \\
        --checkm2_db ${db} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir final_bins/
    touch final_bins_quality_reports.tsv final_bins/bin_1.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """
}