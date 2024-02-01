process BWAMEM2_MEM {
    tag "${meta.id} align to ${ref_fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:6351200
f24497efba12c219c2bea4bb0f69a9d47-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:6351200f24497efba12c219c2bea4bb0f69a9
d47-0' }"


    input:
    tuple val(meta), path(reads), path(ref_fasta), path(ref_fasta_index)

    output:
    tuple val(meta), path("output/${meta.id}_sorted.bam"), path("output/${meta.id}_sorted.bam.bai") , emit: bams
    path  "versions.yml", emit: versions

    script:
    def input_reads = "";
    if ( reads.collect().size() == 1 ) {
        input_reads = "${reads[0]}";
    } else {
        if (reads[0].name.contains("_1")) {
            input_reads = "${reads[0]} ${reads[1]}"
        } else {
            input_reads = "${reads[1]} ${reads[0]}"
        }
    }

    def samtools_args = "-f 12 -F 256 -uS"

    """
    mkdir -p output
    echo "mapping files to reference genome"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_fasta} \
      ${input_reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${meta.id}_sorted.bam

    echo "samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam
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
