process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    val(out_fmt)

    output:
    tuple val(meta), path("*_reformated.${out_fmt}")                                             , emit: reformated
    path  "versions.yml"                                                                         , emit: versions_bbmap_reformat, topic: versions
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions
    path  "*.log"                                                                                , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_reads  = "in=${reads}"
    def out_reads = "out=${prefix}_reformated.${out_fmt}"

    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    reformat.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        threads=${task.cpus} \\
        ${args} \\
        &> ${prefix}.reformat.sh.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_reformated.${out_fmt}
    touch ${prefix}.reformat.sh.log
    """
}
