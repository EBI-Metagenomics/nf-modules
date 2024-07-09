
process KEGGPATHWAYSCOMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kegg-pathways-completeness:1.0.5--pyhdfd78af_0':
        'biocontainers/kegg-pathways-completeness:1.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(filtered_tbl), path(ko_list)

    output:
    tuple val(meta), path("*_contigs.tsv") , emit: kegg_contigs
    tuple val(meta), path("*_pathways.tsv"), emit: kegg_pathways
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kegg_input = ko_list ? "-l ${ko_list}" : "-i ${filtered_tbl}"

    if (ko_list && filtered_tbl){
        log.warn("Both \$ko_list and \$filtered_tbl were given as input types, will fall back to using \$ko_list i.e. ${ko_list}");
    }

    """
    give_pathways \\
    ${kegg_input} \\
    -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_pathways --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.kegg_contigs.tsv
    touch ${prefix}.kegg_pathways.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_pathways --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
