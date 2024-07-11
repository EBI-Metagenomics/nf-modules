process DRAM_DISTILL {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.3.5--pyhdfd78af_0':
        'quay.io/biocontainers/dram:1.3.5--pyhdfd78af_0' }"

    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return "--bind ${moduleDir}/tests/fixtures/dram_distill_dbs/:/data/ --bind ${moduleDir}/tests/fixtures/dram_distill_dbs/CONFIG:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"
        } else {
            return "-v ${moduleDir}/tests/fixtures/dram_distill_dbs/:/data/ -v ${moduleDir}/tests/fixtures/dram_distill_dbs/CONFIG:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"
        }
    }

    input:
    tuple val(meta), path(tsv_input)

    output:
    tuple val(meta), path("*_dram.html"), emit: html
    tuple val(meta), path("*_dram.tsv") , emit: tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.5' // WARN: dram has no option to print the tool version. This is the container version

    """
    DRAM.py \\
        distill \\
        -i ${dram_input}  \\
        -o dram_out \\
        ${args}
    mv dram_out/product.html ${prefix}_dram.html
    mv dram_out/product.tsv ${prefix}_dram.tsv
    
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """
}