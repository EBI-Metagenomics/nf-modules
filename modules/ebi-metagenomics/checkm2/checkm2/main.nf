process CHECKM2 {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bins)
    path checkm_db

    output:
    tuple val(meta), path("*_checkm_output/quality_report.tsv"),  emit: checkm2_stats
    path "versions.yml"                                        ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    checkm2 predict --threads ${task.cpus} \
                    --input ${bins} \
                    --output-directory ${prefix}_checkm_output \
                    --database_path ${checkm_db} \
                    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CheckM2 : \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_checkm_output
    echo "genome\tcompleteness\tcontamination" > ${prefix}_checkm_output/quality_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CheckM2 : \$(checkm2 --version)
    END_VERSIONS
    """
}
