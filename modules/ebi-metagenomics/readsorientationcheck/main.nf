process READSORIENTATIONCHECK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/infernal_seqkit_python:a8730f828bc4151e' :
        'community.wave.seqera.io/library/infernal_seqkit_python:9fddb07b0b9b76de' }"

    input:
    tuple val(meta), path(reads)
    path covariance_model
    val n_reads

    output:
    tuple val(meta), path("*.orientation.tsv"), emit: tsv
    tuple val(meta), path("*.cmsearch.tbl.gz") , emit: cmsearch_tbl
    tuple val(meta), env("classification")    , emit: classification
    tuple val("${task.process}"), val('infernal'), eval("cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//'"), topic: versions, emit: versions_infernal
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^seqkit v//'")                                        , topic: versions, emit: versions_seqkit
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'")                                       , topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Only R1 is screened; for paired-end input R1 is the first file
    def r1 = reads instanceof List ? reads[0] : reads
    """
    # Take the first ${n_reads} R1 reads as FASTA.
    # Number each read so even repeated FASTQ names count independently.
    seqkit head -n ${n_reads} ${r1} \\
        | seqkit fq2fa \\
        | seqkit replace -p '^.+\$' -r '{nr}' \\
        > reads.fa

    if [ -s reads.fa ]; then
        cmsearch \\
            --cpu ${task.cpus} \\
            --noali \\
            --hmmonly \\
            --cut_ga \\
            ${args} \\
            -o /dev/null \\
            --tblout ${prefix}.cmsearch.tbl \\
            ${covariance_model} \\
            reads.fa
    else
        # An empty FASTQ has no hits to parse
        touch ${prefix}.cmsearch.tbl
    fi

    check_r1_orientation.py \\
        --sample ${prefix} \\
        --fasta reads.fa \\
        --tblout ${prefix}.cmsearch.tbl \\
        --cm-database ${covariance_model} \\
        --reads ${n_reads} \\
        ${args2} \\
        > ${prefix}.orientation.tsv

    gzip -n ${prefix}.cmsearch.tbl

    classification=\$(awk -F '\\t' 'NR == 2 { print \$NF }' ${prefix}.orientation.tsv)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.orientation.tsv
    echo "" | gzip > ${prefix}.cmsearch.tbl.gz
    classification="not_mixed"
    """
}
