
process HMMER_HMMSEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(fasta)
    path hmm_db

    output:
    tuple val(meta), path("*.tbl"), emit: tbl
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def hmm_file = "" 
    if (!hmm_db instanceof List) {
        hmm_file = hmm_db
    }
    else{
        for(iter_file in hmm_db) {
            if (iter_file.getExtension() == 'hmm'){
                hmm_file = iter_file
                break
            }
        }
    }

    """
    hmmsearch \\
        ${args} \\
        --cpu ${task.cpus} \\
        --domtblout ${prefix}_hmmsearch.tbl \\
        -o /dev/null \\
        ${hmm_file} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep 'HMMER' | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_hmmsearch.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep 'HMMER' | cut -d' ' -f3)
    END_VERSIONS
    """
}
