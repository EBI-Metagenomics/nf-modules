process SAMTOOLS_BAM2FQ {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
	tuple val(meta), path(bam)
    	tuple val(meta), path(reads)

    output:
    	tuple val(meta), path("*.bwa.fq.gz"), emit: reads
    	path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta[0].id
    def single_end = reads.collect().size() == 1

    if (single_end){
        """
        samtools \\
            bam2fq \\
            -@ $task.cpus \\
            $bam | gzip --no-name > ${prefix}.bwa.fq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        samtools \\
            bam2fq \\
            -@ $task.cpus \\
            -1 ${prefix}_1.bwa.fq.gz \\
            -2 ${prefix}_2.bwa.fq.gz \\
            -0 /dev/null \\
            -s /dev/null \\
            $bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } 
}
