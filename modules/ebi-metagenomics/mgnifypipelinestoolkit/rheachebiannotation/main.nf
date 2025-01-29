process MGNIFYPIPELINESTOOLKIT_RHEACHEBIANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.2.1--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(diamond_tsv)
    path(rhea2chebi)

    output:
    tuple val(meta), path("*.tsv"), emit: rhea2proteins_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_fasta_compressed = fasta.name.endsWith(".gz")
    def is_diamond_compressed = diamond_tsv.name.endsWith(".gz")
    def fasta_file = fasta.name.replace(".gz", "")
    def diamond_file = diamond_tsv.name.replace(".gz", "")
    """
    if [ "$is_fasta_compressed" == "true" ]; then
        gzip -c -d $fasta > ${fasta_file}
    fi

    if [ "$is_diamond_compressed" == "true" ]; then
        gzip -c -d $diamond_tsv > ${diamond_file}
    fi

    add_rhea_chebi_annotation \\
        --diamond_hits ${diamond_file} \\
        --proteins ${fasta_file} \\
        --rhea2chebi ${rhea2chebi} ${args} \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
