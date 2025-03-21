process DRAM_DISTILL {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.3.5--pyhdfd78af_0':
        'biocontainers/dram:1.3.5--pyhdfd78af_0' }"

    containerOptions {
        // If we do anything more fancy to reduce the repeated code the linting fails, the nf-core tools
        // container URL parsing code it's not handling this
        if ( workflow.containerEngine == 'singularity' ) {
            return "--bind ${task.workDir}/${dram_dbs}/:/data/ --bind ${task.workDir}/${dram_dbs}/DRAM_CONFIG.json:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"
        } else {
            return "--volume ${task.workDir}/${dram_dbs}/:/data/ --volume ${task.workDir}/${dram_dbs}/DRAM_CONFIG.json:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"
        }
    }

    input:
    tuple val(meta), path(tsv_input)
    path(dram_dbs)

    output:
    tuple val(meta), path("*_dram.html.gz")               , emit: html
    tuple val(meta), path("*_dram.tsv.gz")                , emit: out_tsv
    tuple val(meta), path("*_genome_stats.tsv.gz")        , emit: stats_tsv
    tuple val(meta), path("*_metabolism_summary.xlsx.gz") , emit: metabolism_xslx
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: dram has no option to print the tool version. This is the container version
    def VERSION = '1.3.5'
    def is_compressed = tsv_input.getExtension() == "gz"
    def tsv_file_name = tsv_input.name.replace(".gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d ${tsv_input} > ${tsv_file_name}
    fi

    DRAM.py \\
        distill \\
        -i ${tsv_input}  \\
        -o dram_out \\
        ${args}

    mv dram_out/product.html ${prefix}_dram.html
    mv dram_out/product.tsv ${prefix}_dram.tsv
    mv dram_out/genome_stats.tsv ${prefix}_genome_stats.tsv
    mv dram_out/metabolism_summary.xlsx ${prefix}_metabolism_summary.xlsx

    gzip ${prefix}_dram.html
    gzip ${prefix}_dram.tsv
    gzip ${prefix}_genome_stats.tsv
    gzip ${prefix}_metabolism_summary.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.5'
    """
    touch ${prefix}_dram.html
    touch ${prefix}_dram.tsv
    touch ${prefix}_genome_stats.tsv
    touch ${prefix}_metabolism_summary.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """
}
