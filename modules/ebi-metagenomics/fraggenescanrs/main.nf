process FRAGGENESCANRS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/frag_gene_scan_rs:1.1.0--h4349ce8_0':
        'biocontainers/frag_gene_scan_rs:1.1.0--h4349ce8_0' }"

    input:
    tuple val(meta), path(fasta)
    val(training_file_name)

    output:
    tuple val(meta), path("*.ffn.gz"), emit: nucleotide_fasta
    tuple val(meta), path("*.faa.gz"), emit: amino_acid_fasta
    tuple val(meta), path("*.out.gz"), emit: gene_annotations
    tuple val(meta), path("*.gff.gz"), emit: gff
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def catcommand = fasta.name.endsWith(".gz") ? "zcat" : "cat"
    """
    ${catcommand} ${fasta} | FragGeneScanRs \\
        $args \\
        --training-file ${training_file_name} \\
        --thread-num ${task.cpus} \\
        --output-prefix ${prefix}

    gzip *.{ffn,faa,out,gff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FragGeneScanRs: \$(FragGeneScanRs --version | sed 's/FragGeneScanRs //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.out
    touch ${prefix}.gff

    gzip *.{ffn,faa,out,gff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FragGeneScanRs: \$(FragGeneScanRs --version | sed 's/FragGeneScanRs //')
    END_VERSIONS
    """
}
