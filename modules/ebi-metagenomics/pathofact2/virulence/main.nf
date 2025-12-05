process PATHOFACT2_VIRULENCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/microbiome-informatics/pathofact2_env:v1.0.0'

    input:
    tuple val(meta), path(faa)
    path pathofact_db

    output:
    tuple val(meta), path("*_classifier_virulence_filtered.tsv"), emit: tsv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vf_predict.py \\
        --file ${faa} \
        --model ${pathofact_db}/VF/final_model.joblib \
        --cpus ${task.cpus} \
        --outfile ${prefix}_classifier_virulence.tsv 

    awk -v ${args} -F'\\t' 'NR==1 || \$3 > thr { print }' ${prefix}_classifier_virulence.tsv > ${prefix}_classifier_virulence_filtered.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}_classifier_virulence_filtered.tsv
    """
}
