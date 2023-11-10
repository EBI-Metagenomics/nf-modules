process ESMFOLD {
    tag "$meta.id"
    label 'process_high'

    conda '/home/vangelis/miniconda3/envs/esmfold_new_test' // /hps/nobackup/rdf/metagenomics/service-team/users/vangelis/miniconda3/envs/esmfold_new

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*pdb"), emit: pdb
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    esm-fold \\
        -i ${fasta_name} \\
        -o ${prefix} \\
        --cpu-only \\
        ${args}

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ESM-2: v1.0.3
    END_VERSIONS
    """
}
