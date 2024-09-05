process CRISPRCASFINDER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'microbiome-informatics/genomes-pipeline.crisprcasfinder:4.3.2':
        'microbiome-informatics/genomes-pipeline.crisprcasfinder:4.3.2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("crisprcasfinder_results/GFF/*.gff"), emit: gff
    tuple val(meta), path("crisprcasfinder_results/TSV"), emit: tsv
    tuple val(meta), path("crisprcasfinder_results/rawCRISPRs.fna"), optional: true, emit: fna
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    CRISPRCasFinder.pl -i $fasta \
        ${args} \\
        -so /opt/CRISPRCasFinder/sel392v2.so \
        -outdir crisprcasfinder_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CRISPRCasFinder: \$(echo \$(CRISPRCasFinder.pl -v) | grep -o "version [0-9.]*" | sed "s/version //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p crisprcasfinder_results/TSV
    mkdir -p crisprcasfinder_results/GFF
    touch crisprcasfinder_results/GFF/${prefix}_crisprcasfinder.gff
    touch crisprcasfinder_results/rawCRISPRs.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CRISPRCasFinder: \$(echo \$(CRISPRCasFinder.pl -v) | grep -o "version [0-9.]*" | sed "s/version //g")
    END_VERSIONS
    """
}
