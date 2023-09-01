process FETCHTOOL_ASSEMBLY {
    tag "$assembly_accession"

    label 'process_single'

    container "microbiome-informatics/fetch-tool:v0.9.0"

    input:
    tuple val(meta), val(assembly_accession)
    path fetchtool_config

    output:
    tuple val(meta), path("download_folder/*/raw/${assembly_accession}*.fasta.gz"), emit: assembly
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$assembly_accession"
    """
    fetch-assembly-tool -d download_folder/ -as $assembly_accession -c $fetchtool_config -v $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch-tool: \$(fetch-assembly-tool --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$assembly_accession"
    """
    mkdir -p download_folder/test/raw/

    touch download_folder/test/raw/${assembly_accession}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch-tool: \$(fetch-assembly-tool --version)
    END_VERSIONS
    """
}
