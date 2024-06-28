
process ANTISMASH {
    tag "$meta.id"
    label 'process_medium'

    // current antiSMASH bioconda container is buggy so we use our custom one (https://github.com/antismash/antismash/discussions/711)
    container "microbiome-informatics/antismash:7.1.0.1_2"

    input:
    tuple val(meta)         , path(contigs)
    tuple val(meta)         , path(genes)
    tuple val(meta)         , path(gbk_input)
    tuple path(antismash_db), val(db_version)

    output:
    tuple val(meta), path("${meta.id}_results/${meta.id}.gbk")  , emit: gbk
    tuple val(meta), path("${meta.id}_results/${meta.id}.json") , emit: json
    tuple val(meta), path("${meta.id}_antismash.tar.gz")        , emit: results_tarball
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (gbk_input == []){
        """
        antismash \\
            -c ${task.cpus} \\
            ${args} \\
            --databases ${antismash_db} \\
            --genefinding-gff3 ${genes} \\
            --output-dir ${prefix}_results \\
            ${contigs}

        tar -czf ${prefix}_antismash.tar.gz ${prefix}_results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
            antiSMASH database: ${db_version}
        END_VERSIONS
        """
    }
    else{
        """
        antismash \\
            -c ${task.cpus} \\
            ${args} \\
            --databases ${antismash_db} \\
            --output-dir ${prefix}_results \\
            ${gbk_input}

        tar -czf ${prefix}_antismash.tar.gz ${prefix}_results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
            antiSMASH database: ${db_version}
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_antismash.tar.gz
    mkdir ${prefix}_results/
    touch ${prefix}_results/${prefix}.gbk
    touch ${prefix}_results/${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
        antiSMASH database: ${db_version}
    END_VERSIONS
    """
}
