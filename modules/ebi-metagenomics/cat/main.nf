process CAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cat:5.3--hdfd78af_0':
        'biocontainers/cat:5.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(contigs_fasta)
    path cat_db
    path taxonomy_db

    output:
    tuple val(meta), path('*_summary.txt'), emit: summary
    tuple val(meta), path('*_contig2classification_official_names.tsv'), emit: contigs_classification
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Analysing contigs"
    CAT contigs \
        -n ${task.cpus} \
        -c ${contigs_fasta} \
        -d ${cat_db} \
        -t ${taxonomy_db} \
        --out_prefix ${prefix}

    echo "Adding taxonomy names"
    CAT add_names \
        -i ${prefix}.contig2classification.txt \
        -o ${prefix}_contig2classification_official_names.tsv \
        -t ${taxonomy_db} --only_official

    echo "Summarizing output"
    CAT summarise -c ${contigs_fasta} \
        -i ${prefix}_contig2classification_official_names.tsv \
        -o ${prefix}_summary.txt

    ############################################################
    # Storage optimization                                     #
    # Getting rid of files that are not needed after it runs   #
    ############################################################
    rm -f "${prefix}.alignment.diamond"
    rm -f "${prefix}.predicted_proteins.faa"
    rm -f "${prefix}.predicted_proteins.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_contig2classification_official_names.tsv ${prefix}_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
    END_VERSIONS
    """
}
