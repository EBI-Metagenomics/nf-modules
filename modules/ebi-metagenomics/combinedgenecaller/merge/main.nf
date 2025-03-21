process COMBINEDGENECALLER_MERGE {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.0.1--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pyrodigal_gff, stageAs: "pyrodigal/"), path(pyrodigal_ffn, stageAs: "pyrodigal/"), path(pyrodigal_faa, stageAs: "pyrodigal/"), path(fgs_gff, stageAs: "fgsrs/"), path(fgs_ffn, stageAs: "fgsrs/"), path(fgs_faa, stageAs: "fgsrs/"), path(mask)

    output:
    tuple val(meta), path("*.faa.gz"), emit: faa
    tuple val(meta), path("*.ffn.gz"), emit: ffn
    tuple val(meta), path("*.gff.gz"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Check if files are compressed
    def is_prod_gff_compressed = pyrodigal_gff.name.endsWith(".gz")
    def is_prod_ffn_compressed = pyrodigal_ffn.name.endsWith(".gz")
    def is_prod_faa_compressed = pyrodigal_faa.name.endsWith(".gz")
    def is_fgs_gff_compressed = fgs_gff.name.endsWith(".gz")
    def is_fgs_ffn_compressed = fgs_ffn.name.endsWith(".gz")
    def is_fgs_faa_compressed = fgs_faa.name.endsWith(".gz")

    // Get uncompressed file names
    def pyrodigal_gff_file = pyrodigal_gff.name.replace(".gz", "")
    def pyrodigal_ffn_file = pyrodigal_ffn.name.replace(".gz", "")
    def pyrodigal_faa_file = pyrodigal_faa.name.replace(".gz", "")
    def fgs_gff_file = fgs_gff.name.replace(".gz", "")
    def fgs_ffn_file = fgs_ffn.name.replace(".gz", "")
    def fgs_faa_file = fgs_faa.name.replace(".gz", "")

    def is_mask_compressed = false
    def mask_parameter = ""
    if ( mask && mask.size() > 0 ) {
        is_mask_compressed = mask.name.endsWith(".gz")
        mask_parameter = "--mask ${mask.name.replace(".gz", "")} "
    }
    """
    if [ "${is_prod_gff_compressed}" == "true" ]; then
        gzip -d -c ${pyrodigal_gff} > ${pyrodigal_gff_file}
    fi
    if [ "${is_prod_ffn_compressed}" == "true" ]; then
        gzip -d -c ${pyrodigal_ffn} > ${pyrodigal_ffn_file}
    fi
    if [ "${is_prod_faa_compressed}" == "true" ]; then
        gzip -d -c ${pyrodigal_faa} > ${pyrodigal_faa_file}
    fi
    if [ "${is_fgs_gff_compressed}" == "true" ]; then
        gzip -d -c ${fgs_gff} > ${fgs_gff_file}
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

    combined_gene_caller_merge \\
        ${args} \\
        -n ${prefix} \\
        --pyrodigal-gff ${pyrodigal_gff_file} \\
        --pyrodigal-ffn ${pyrodigal_ffn_file} \\
        --pyrodigal-faa ${pyrodigal_faa_file} \\
        --fgsrs-gff ${fgs_gff_file} \\
        --fgsrs-ffn ${fgs_ffn_file} \\
        --fgsrs-faa ${fgs_faa_file} ${mask_parameter}

    gzip -n ${prefix}.{faa,ffn,gff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.{faa,ffn,gff}

    gzip -n ${prefix}.{faa,ffn,gff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
