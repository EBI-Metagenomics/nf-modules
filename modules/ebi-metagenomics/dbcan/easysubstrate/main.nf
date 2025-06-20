process DBCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0'
        : 'biocontainers/dbcan:5.1.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gff) // only used for cluster with a proteins file
    tuple path(dbcan_db, stageAs: "dbcan_db"), val(db_version)
    val(mode)

    output:
    tuple val(meta), path("results/${prefix}_overview.tsv.gz"),                             emit: overview_txt
    tuple val(meta), path("results/${prefix}_dbcan_hmm_results.tsv.gz"),                    emit: dbhmm_output_tsv
    tuple val(meta), path("results/${prefix}_dbcansub_hmm_results.tsv.gz"),                 emit: dbsub_output_tsv
    tuple val(meta), path("results/${prefix}_diamond.out.gz"),                              emit: diamond_output
    tuple val(meta), path("results/${prefix}_uniinput.faa.gz"),             optional: true, emit: uniinput_faa
    tuple val(meta), path("results/${prefix}_cgc.gff.gz"),                  optional: true, emit: cgc_gff
    tuple val(meta), path("results/${prefix}_cgc_standard_out.tsv.gz"),     optional: true, emit: cgc_standard_tsv
    tuple val(meta), path("results/${prefix}_substrate_prediction.tsv.gz"), optional: true, emit: substrate_prediction_tsv
    tuple val(meta), path("results/synteny_pdf/*-syntenic.pdf.gz"),         optional: true, emit: synteny_pdfs
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def is_fasta_compressed = fasta.getExtension() == "gz"
    def fasta_name = is_fasta_compressed ? fasta.getBaseName() : fasta

    def is_gff_compressed = gff.getExtension() == "gz"
    def gff_name = is_gff_compressed ? gff.getBaseName() : gff
    """
    if [ "${is_fasta_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if [ "${is_gff_compressed}" == "true" ]; then
        gzip -c -d ${gff} > ${gff_name}
    fi

    run_dbcan \\
        easy_substrate \\
        ${args} \\
        --threads ${task.cpus} \\
        --db_dir dbcan_db \\
        --output_dir results \\
        --input_raw_data ${fasta_name} \\
        --gff_type prodigal \\
        --input_gff ${gff_name} \\
        --mode ${mode}

    # Bulk rename of the results, dbcan doesn't have a prefix parameter
    find results -type f | while read -r file; do
        mv "\$file" "\$(dirname "\$file")/${prefix}_\$(basename "\$file" | tr [:upper:] [:lower:])"
    done

    gzip results/*.*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(run_dbcan version | sed "s/dbCAN version: //")
        dbcan_db: "${db_version}"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir results
    touch results/${prefix}_overview.tsv
    touch results/${prefix}_dbcan_hmm_results.tsv
    touch results/${prefix}_dbCANsub_hmm_results.tsv
    touch results/${prefix}_diamond.out

    gzip results/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(run_dbcan version | sed "s/dbCAN version: //")
        dbcan_db: ${db_version}
    END_VERSIONS
    """
}
