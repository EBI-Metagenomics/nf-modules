process ESMFOLD {
    tag "$meta.id"
    label 'process_high'

    conda params.esm_conda_path

    input:
    tuple val(meta), path(fasta)
    val(compute_mode)

    output:
    tuple val(meta), path("${meta.id}")           , emit: pdb
    tuple val(meta), path("${meta.id}_scores.txt"), emit: scores
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def compute_flag = compute_mode == 'gpu' ? '' : '--cpu-only'
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    esm-fold \\
        -i ${fasta_name} \\
        -o ${prefix} \\
        ${compute_flag} \\
        ${args} > ${prefix}_scores.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ESM-2: v1.0.3
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${fasta}.pdb
    touch ${prefix}_scores.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ESM-2: v1.0.3
    END_VERSIONS
    """
}
