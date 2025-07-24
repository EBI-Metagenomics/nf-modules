process BINNET {
    tag "$meta ? $meta.id : $meta2.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ 'community.wave.seqera.io/library/binette:1.1.2--90fa8b078479d111' }"

    input:
    tuple val(meta), path(fasta_files, stageAs: "input_bins/*")
    tuple val(meta2), path(tsv_files, stageAs: "input_tsv/*")
    tuple val(meta3), path(contigs)
    tuple val(meta4), path(proteins)
    path db

    output:
    tuple val(meta), path("final_bins/*")                  , emit: refined_bins
    tuple val(meta), path("final_bins_quality_reports.tsv"), emit: refined_bins_report
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args              = task.ext.args ?: ''
        bin_dirs          = fasta_files ? "--bin_dirs input_bins/*" : ""
        contig2bin_tables = tsv_files & !fasta_files ? "--contig2bin_tables input_tsv/*" : ""
        proteins          = proteins ? "--proteins ${proteins}" : ""
    """
    binnete \\
        ${bin_dirs} \\
        ${contig2bin_tables} \\
        ${proteins} \\
        --contigs ${contigs} \\
        --output . \\
        --threads ${task.cpus} \\
        --checkm2_db ${db} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binnete: \$( binnete --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir final_bins/
    touch final_bins_quality_reports.tsv final_bins/bin_1.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binnete: \$( binnete --version )
    END_VERSIONS
    """
}