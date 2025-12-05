process PATHOFACT2_TOXINS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/microbiome-informatics/pathofact2_env:v1.0.0'

    input:
    tuple val(meta), path(faa)
    path pathofact_db

    output:
    tuple val(meta), path("*_classifier_toxins_filtered.tsv"), emit: tsv
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // This customizes the command: toxins prediction
    def args2 = task.ext.args2 ?: ''
    // This customizes the command: awk filtering
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tox_predict.py \\
        ${args} \\
        -s ${faa} \\
        -m ${pathofact_db}/TOX/final_model.joblib \\
        -v ${model_tox}/TOX/ \\
        --cpus ${task.cpus} \\
        -o ${prefix}_classifier_toxins.tsv
 
    awk -v ${args2} -F'\\t' 'NR==1 || \$3 > thr { print }' ${prefix}_classifier_toxins.tsv > ${prefix}_classifier_toxins_filtered.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    echo $args2
    touch ${prefix}_classifier_toxins_filtered.tsv
    """
}
