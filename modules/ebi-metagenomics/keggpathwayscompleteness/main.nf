
process KEGGPATHWAYSCOMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kegg-pathways-completeness:1.0.3--pyhdfd78af_0':
        'biocontainers/kegg-pathways-completeness:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(hmmsearch_tbl)
    tuple val(meta), path(ko_list)

    output:
    tuple val(meta), path("*_contigs.tsv") , emit: kegg_contigs
    tuple val(meta), path("*_pathways.tsv"), emit: kegg_pathways
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = "1.0.3" // No way to get version automatically so hard-coding it

    if (ko_list == []){
        """
        sed \\
        '/^#/d; s/ \\+/\\t/g' \\
        ${hmmsearch_tbl} \\
        > ${prefix}.filtered.tbl

        parsing_hmmscan \\
        -i ${prefix}.filtered.tbl \\
        -f ${fasta}

        give_pathways \\
        -i ${prefix}.filtered.tbl_parsed \\
        -o ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            keggpathwayscompleteness: ${version}
        END_VERSIONS
        """
    }
    else {
        """
        give_pathways \\
        -l ${ko_list} \\
        -o ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            keggpathwayscompleteness: ${version}
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = "1.0.3" // No way to get version automatically so hard-coding it

    """
    touch ${prefix}.kegg_contigs.tsv
    touch ${prefix}.kegg_pathways.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        keggpathwayscompleteness: ${version}
    END_VERSIONS
    """
}
