process BWAMEM2_MEM {
    label 'process_high'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(ref_index)
    val align // if true: align (include reads), else: decontaminate (exclude reads)


    output:
    tuple val(meta), path("*_sorted.bam"), path("*_sorted.bam.bai"), emit: bam
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = "-M"
    def prefix = task.ext.prefix ?: meta[0].id
    def database = task.ext.database ?: meta2[0].id

    def samtools_args = ""
    if ( align ) {
        args2 = "-q 20 -Sb"
    } else {
        args2 = "-f 12 -F 256 -uS"
    }

    """
    bwa-mem2 mem \\
        $args \\
        -t $task.cpus \\
        $database \\
        $reads \\
        | samtools view -@ ${task.cpus} $args2 - \\
        | samtools sort -@ ${task.cpus} -O bam - -o ${prefix}_sorted.bam
    samtools index -@ ${task.cpus} ${prefix}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
