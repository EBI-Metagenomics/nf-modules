process PATHOFACT2_TOXINS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "microbiome-informatics/pathofact2_env:v1.0.2"

    input:
    tuple val(meta), path(fasta)
    path zenodo_file

    output:
    tuple val(meta), path("*_classifier_toxins_filtered.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''               // This customizes the command: toxins prediction
    def args2 = task.ext.args2 ?: 'thr=0.6'      // This customizes the command: awk filtering
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''
    """
    $uncompress_input

    tar -xavf ${zenodo_file}

    tox_predict.py \\
        ${args} \\
        -s ${fasta_name} \\
        -m Models/TOX/final_model.joblib \\
        -v Models/TOX/ \\
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
