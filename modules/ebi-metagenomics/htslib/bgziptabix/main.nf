process HTSLIB_BGZIPTABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.24--ha79157c_0' :
        'quay.io/biocontainers/htslib:1.24--ha79157c_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz")    , emit: gz
    tuple val(meta), path("*.gz.gzi"), emit: gzi
    tuple val(meta), path("*.gz.csi"), emit: csi
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bgzip
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | sed -n 's/^Version: *//p'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_gz = input.getExtension() == "gz"
    // keep the format extension (tabix uses it to detect the preset), e.g. "x.gff.gz" -> "gff", "x.gff" -> "gff"
    def extension = in_gz ? input.getBaseName().tokenize(".")[-1] : input.getExtension()
    def output = "${prefix}.${extension}.gz"
    // a staged input with the same name as the output would be overwritten
    if ("${input}" == output) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    // a compressed input is decompressed first, because it can be gzip and not BGZF
    def reader = in_gz ? "bgzip --threads ${task.cpus} --decompress --stdout ${input}" : "cat ${input}"
    """
    # bgzip with a .gzi index, for paginated access to the compressed file
    ${reader} | bgzip \\
        --threads ${task.cpus} \\
        --stdout \\
        --index \\
        --index-name ${output}.gzi \\
        ${args} \\
        > ${output}

    # tabix with a .csi index, for access by genomic region (it scales better than .tbi for large contigs)
    tabix \\
        --threads ${task.cpus} \\
        --csi \\
        ${args2} \\
        ${output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_gz = input.getExtension() == "gz"
    def extension = in_gz ? input.getBaseName().tokenize(".")[-1] : input.getExtension()
    def output = "${prefix}.${extension}.gz"
    if ("${input}" == output) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo "" | gzip > ${output}
    touch ${output}.gzi
    touch ${output}.csi
    """
}
