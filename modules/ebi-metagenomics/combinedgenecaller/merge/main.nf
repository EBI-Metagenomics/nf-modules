process COMBINEDGENECALLER_MERGE {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.2.1--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(prodigal_sco, stageAs: "prodigal/"), path(prodigal_ffn, stageAs: "prodigal/"), path(prodigal_faa, stageAs: "prodigal/"), path(fgs_out, stageAs: "fgsrs/"), path(fgs_ffn, stageAs: "fgsrs/"), path(fgs_faa, stageAs: "fgsrs/"), path(mask)

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

    // Check if files are compressed
    def is_prod_sco_compressed = prodigal_sco.name.endsWith(".gz")
    def is_prod_ffn_compressed = prodigal_ffn.name.endsWith(".gz")
    def is_prod_faa_compressed = prodigal_faa.name.endsWith(".gz")
    def is_fgs_out_compressed = fgs_out.name.endsWith(".gz")
    def is_fgs_ffn_compressed = fgs_ffn.name.endsWith(".gz")
    def is_fgs_faa_compressed = fgs_faa.name.endsWith(".gz")

    // Get uncompressed file names
    def prodigal_sco_file = prodigal_sco.name.replace(".gz", "")
    def prodigal_ffn_file = prodigal_ffn.name.replace(".gz", "")
    def prodigal_faa_file = prodigal_faa.name.replace(".gz", "")
    def fgs_out_file = fgs_out.name.replace(".gz", "")
    def fgs_ffn_file = fgs_ffn.name.replace(".gz", "")
    def fgs_faa_file = fgs_faa.name.replace(".gz", "")

    def is_mask_compressed = false
    def mask_parameter = ""
    if ( mask ) {
        is_mask_compressed = mask.name.endsWith(".gz")
        mask_parameter = "--mask ${mask.name.replace(".gz", "")} "
    }
    """
    if [ "${is_prod_sco_compressed}" == "true" ]; then
        gzip -d -c ${prodigal_sco} > ${prodigal_sco_file}
    fi
    if [ "${is_prod_ffn_compressed}" == "true" ]; then
        gzip -d -c ${prodigal_ffn} > ${prodigal_ffn_file}
    fi
    if [ "${is_prod_faa_compressed}" == "true" ]; then
        gzip -d -c ${prodigal_faa} > ${prodigal_faa_file}
    fi
    if [ "${is_fgs_out_compressed}" == "true" ]; then
        gzip -d -c ${fgs_out} > ${fgs_out_file}
    fi
    if [ "${is_fgs_ffn_compressed}" == "true" ]; then
        gzip -d -c ${fgs_ffn} > ${fgs_ffn_file}
    fi
    if [ "${is_fgs_faa_compressed}" == "true" ]; then
        gzip -d -c ${fgs_faa} > ${fgs_faa_file}
    fi
    if [[ "${is_mask_compressed}" == "true" ]]; then
        gunzip ${mask}
    fi

    cgc_merge \\
        ${args} \\
        -n ${prefix} \\
        --prodigal-out ${prodigal_sco_file} \\
        --prodigal-ffn ${prodigal_ffn_file} \\
        --prodigal-faa ${prodigal_faa_file} \\
        --fgs-out ${fgs_out_file} \\
        --fgs-ffn ${fgs_ffn_file} \\
        --fgs-faa ${fgs_faa_file} ${mask_parameter}

    gzip -n ${prefix}.{faa,ffn,out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.{faa,ffn,out}

    gzip -n ${prefix}.{faa,ffn,out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
