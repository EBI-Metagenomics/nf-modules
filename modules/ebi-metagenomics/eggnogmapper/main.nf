process EGGNOGMAPPER {
    tag "$meta.id"
    label 'process_long'

    conda "bioconda::eggnog-mapper=2.1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(annotation_hit_table)
    path(eggnog_data_dir)
    path(eggnog_db)
    path(eggnog_diamond_db)

    output:
    tuple val(meta), path("*.hits")           , emit: hits, optional: true
    tuple val(meta), path("*.annotations")    , emit: annotations, optional: true
    tuple val(meta2), path("*.seed_orthologs"), emit: orthologs, optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def annotation_hit_input = annotation_hit_table ? "--annotate_hits_table ${annotation_hit_table}" : ""
    def eggnog_db_input = eggnog_db ? "--database ${eggnog_db}" : ""
    def eggnog_diamond_db_input = eggnog_diamond_db ? "--dmnd_db ${eggnog_diamond_db}" : ""
    def dbmem = task.memory.toMega() > 40000 ? '--dbmem' : ''

    def fasta_bool = fasta ? fasta.name : "no_fasta"
    def is_compressed = fasta_bool.endsWith(".gz")
    def fasta_name = fasta_bool.replace(".gz", "")
    def fasta_input = fasta ? "-i ${fasta_name}" : ""

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    emapper.py \\
        ${args} \\
        ${fasta_input} \\
        ${annotation_hit_input} \\
        ${eggnog_db_input} \\
        ${eggnog_diamond_db_input} \\
        ${dbmem} \\
        --cpu ${task.cpus} \\
        --data_dir ${eggnog_data_dir} \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emapper.hits

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
    END_VERSIONS
    """
}
