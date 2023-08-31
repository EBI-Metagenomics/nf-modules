process ALIGNMENT {

    tag "${meta.id} align to ${ref_fasta}"

    label 'alignment'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(input_ch_reads)
    tuple val(meta), path(ref_fasta), path(ref_fasta_index)

    output:
    tuple val(meta), path(ref_fasta), path("output/${meta.id}_sorted.bam"), path("output/${meta.id}_sorted.bam.bai") , emit: bams

    script:

    // define reads
    reads = input_ch_reads.collect()
    def input_reads = "";
    if ( meta.single_end ) {
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

    echo "mapping files to host genome"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_fasta} \
      ${input_reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${meta.id}_sorted.bam

    echo "samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam
    """
}
