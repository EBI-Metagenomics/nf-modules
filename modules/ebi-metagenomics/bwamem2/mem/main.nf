process BWAMEM2_MEM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:6351200f24497efba12c219c2bea4bb0f69a9d47-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:6351200f24497efba12c219c2bea4bb0f69a9d47-0' }"
    

    input:
        tuple val(meta), path(reads)
        path(index_dir)

    output:
        tuple val(meta), path("${meta.id}_sorted.bam")    , emit: bam 
        tuple val(meta), path("${meta.id}_sorted.bam.bai"), emit: bai
        path  "versions.yml"                              , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = "-M"
    def args2 = "-f 12 -F 256 -uS"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    script:
    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools view -@ ${task.cpus} $args2 - \\
        | samtools sort -@ ${task.cpus} -O bam - -o ${prefix}_sorted.bam
    samtools index -@ ${task.cpus} ${prefix}_sorted.bam
    """

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
