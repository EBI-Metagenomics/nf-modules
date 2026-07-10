
process GENOMEPROPERTIES {
    tag "$meta.id"
    label 'process_single'

    container 'microbiome-informatics/genome-properties:v2.0.2'

    input:
    tuple val(meta), path(ips)

    output:
    tuple val(meta), path("*.txt.gz") , emit: summary
    tuple val(meta), path("*.json.gz"), emit: json
    tuple val(meta), path("*.tsv.gz") , emit: tsv
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gp_version = "2.0.2" // No way to get the version from the tool directly so have to hardcode

    """
    assign_genome_properties.pl \\
        ${args} \\
        -matches ${ips} \\
        -gpdir /opt/genome-properties/flatfiles/ \\
        -gpff genomeProperties.txt \\
        -name ${prefix}

    mv JSON_${prefix} ${prefix}_gp.json
    mv SUMMARY_FILE_${prefix} ${prefix}_gp.txt
    mv TABLE_${prefix} ${prefix}_gp.tsv

    gzip ${prefix}_gp.{json,txt,tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Genome Properties: ${gp_version}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gp_version = "2.0.2" // No way to get the version from the tool directly so have to hardcode

    """
    touch ${prefix}_gp.json
    touch ${prefix}_gp.txt
    touch ${prefix}_gp.tsv

    gzip ${prefix}_gp.{json,txt,tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Genome Properties: ${gp_version}
    END_VERSIONS
    """
}
