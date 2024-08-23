
process FASTQSUFFIXHEADERCHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.5--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end){
        """
        fastq_suffix_header_check -f 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqsuffixheadercheck: \$(samtools --version |& sed '1!d ; s/samtools //')
        END_VERSIONS
        """

    } else{
        """
        samtools \\
            sort \\
            $args \\
            -@ $task.cpus \\
            -o ${prefix}.bam \\
            -T $prefix \\
            $bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqsuffixheadercheck: \$(samtools --version |& sed '1!d ; s/samtools //')
        END_VERSIONS
        """

    }



    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqsuffixheadercheck: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
