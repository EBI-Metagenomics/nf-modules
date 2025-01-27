process COMBINEDGENECALLER_MERGE {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.2.0--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(prodigal_sco, stageAs: "pro.sco"), path(prodigal_ffn, stageAs: "pro.ffn"), path(prodigal_faa, stageAs: "pro.faa"), path(fgs_out, stageAs: "fgf.out"), path(fgs_ffn, stageAs: "fgf.ffn"), path(fgs_faa, stageAs: "fgs.faa"), path(mask_file)

    output:
    tuple val(meta), path("*.faa.gz"), emit: faa
    tuple val(meta), path("*.ffn.gz"), emit: ffn
    tuple val(meta), path("*.out.gz"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mask_param = mask_file ? " --mask mask_file_uncompressed.txt " : ""
    """
    gzip -cdf ${prodigal_sco} > prodigal_uncompressed.sco
    gzip -cdf ${prodigal_ffn} > prodigal_uncompressed.ffn
    gzip -cdf ${prodigal_faa} > prodigal_uncompressed.faa
    gzip -cdf ${fgs_out} > fgs_uncompressed.out
    gzip -cdf ${fgs_ffn} > fgs_uncompressed.ffn
    gzip -cdf ${fgs_faa} > fgs_uncompressed.faa
    if [[ -n "${mask_file}" ]]; then
        gzip -cdf ${mask_file} > mask_file_uncompressed.txt
    fi

    cgc_merge \\
        ${args} \\
        -n ${prefix} \\
        --prodigal-out prodigal_uncompressed.sco \\
        --prodigal-ffn prodigal_uncompressed.ffn \\
        --prodigal-faa prodigal_uncompressed.faa \\
        --fgs-out fgs_uncompressed.out \\
        --fgs-ffn fgs_uncompressed.ffn \\
        --fgs-faa fgs_uncompressed.faa ${mask_param}

    gzip -n ${prefix}.{faa,ffn,out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa
    touch ${prefix}.ffn

    gzip ${prefix}.faa
    gzip ${prefix}.ffn

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
