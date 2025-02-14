process SANNTIS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'biocontainers/sanntis:0.9.4.0--pyhdfd78af_0'
        : 'biocontainers/sanntis:0.9.4.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(interproscan), path(gbk), path(faa)

    output:
    tuple val(meta), path("*_sanntis.gff.gz"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (faa && gbk) {
        error("SanntiS supports either a GBK or a FAA as the secondary input.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_ips_compressed = interproscan.name.endsWith(".gz")
    def interproscan_file = interproscan ? interproscan.name.replace(".gz", "") : is_ips_compressed

    // Handle the GBK or FAA
    def is_gbk = false
    def input_file = ""
    def decompress = ""
    if (gbk) {
        is_gbk = true
        input_file = gbk.name
        if (gbk.name.endsWith(".gz")) {
            input_file = input_file.replace(".gz", "")
            decompress = "gzip -c -d ${gbk} > ${input_file}"
        }
    }
    if (faa) {
        input_file = faa.name
        if (faa.name.endsWith(".gz")) {
            input_file = input_file.replace(".gz", "")
            decompress = "gzip -c -d ${faa} > ${input_file}"
        }
    }
    def input_arg = is_gbk ? "${prefix}_prepped.gbk" : "--is_protein ${input_file}"
    """
    if [ "${is_ips_compressed}" == "true" ]; then
         gzip -c -d ${interproscan} > ${interproscan_file}
    fi

    $decompress

    if [ "${is_gbk}" == "true" ]; then
        grep -v "/protein_id=" ${input_file} > ${prefix}_prepped.gbk
    fi

    sanntis \\
        --ip-file ${interproscan_file} \\
        --outfile ${prefix}_sanntis.gff \\
        --cpu ${task.cpus} \\
        ${args} \\
        ${input_arg}

    gzip ${prefix}_sanntis.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sanntis.gff
    gzip ${prefix}_sanntis.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
    END_VERSIONS
    """
}
