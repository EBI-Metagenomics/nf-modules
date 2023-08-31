
process FRAGGENESCAN {

    tag "$meta.id"

    label 'process_low'

    conda "bioconda::fraggenescan=1.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fraggenescan:1.31--hec16e2b_4':
        'biocontainers/fraggenescan:1.31--hec16e2b_4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.ffn.gz"), emit: nucleotide_fasta
    tuple val(meta), path("*.faa.gz"), emit: amino_acid_fasta
    tuple val(meta), path("*.out.gz"), emit: gene_annotations
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    args += args.contains("-complete") ? "" : " -complete=0 "
    args += args.contains("-train") ? "" : " -train=illumina_10 "

    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")

    def VERSION = '1.31' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    run_FragGeneScan.pl \
        -thread=$task.cpus \
        -out=$prefix \
        -genome=$fasta_name \
        $args

    gzip *.{ffn,faa,out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fraggenescan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '1.31' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fraggenescan: $VERSION
    END_VERSIONS
    """
}
