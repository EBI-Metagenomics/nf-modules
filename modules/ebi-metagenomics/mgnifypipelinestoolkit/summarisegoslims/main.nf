
process MGNIFYPIPELINESTOOLKIT_SUMMARISEGOSLIMS {
    tag "$meta.id"
    label 'process_single'

    container 'microbiome-informatics/mgnify-pipelines-toolkit:1.0.4'

    input:
    tuple val(meta), path(interproscan_tsv)
    tuple val(meta2), path(gaf)
    path go_obo
    path go_banding

    output:
    tuple val(meta), path("*_summary.tsv"),      emit: go_summary
    tuple val(meta), path("*_summary_slim.tsv"), emit: goslim_summary
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def interproscan_tsv_file = interproscan_tsv.name.replace(".gz", "")
    """
    if [ "${interproscan_tsv.name.endsWith(".gz")}" == "true" ]; then
        gzip -c -d $interproscan_tsv > $interproscan_tsv_file
    fi
    summarise_goslims \\
        ${args} \\
        -go ${go_obo} \\
        -gb ${go_banding} \\
        -i ${interproscan_tsv_file} \\
        -gaf ${gaf} \\
        -o ${prefix}_summary

    mv ${prefix}_summary ${prefix}_summary.tsv
    mv ${prefix}_summary_slim ${prefix}_summary_slim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.tsv
    touch ${prefix}_summary_slim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
