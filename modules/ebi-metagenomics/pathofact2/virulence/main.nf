process PATHOFACT2_VIRULENCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "microbiome-informatics/pathofact2_env:v1.0.4"

    input:
    tuple val(meta), path(fasta)
    path zenodo_file

    output:
    tuple val(meta), path("*_classifier_virulence_filtered.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'thr=0.9'     // This customizes the command: awk filtering
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''
    """
    $uncompress_input

    tar -xavf ${zenodo_file}

    vf_prediction2.py \\
        --file ${fasta_name} \
        --model Models/VF/final_model.joblib \
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
