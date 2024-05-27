
process ANTISMASH {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash:7.1.0.1_2':
        'quay.io/microbiome-informatics/antismash:7.1.0.1_2' }"

    container 'quay.io/microbiome-informatics/antismash:7.1.0.1_2'

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(genes)
    tuple path(antismash_db), val(db_version)

    output:
    tuple val(meta), path("${meta.id}_results/${meta.id}.gbk")  , emit: gbk
    tuple val(meta), path("${meta.id}_antismash.tar.gz")        , emit: results_tarball
    tuple val(meta), path("${meta.id}_regions.json")            , emit: regions
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = antismash_db.name.endsWith(".tar.gz")
    def antismashdb_name = antismash_db.name.replace(".tar.gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        tar -xzvf $antismash_db
    fi

    antismash \\
        -c ${task.cpus} \\
        ${args} \\
        --databases ${antismashdb_name} \\
        --genefinding-gff3 ${genes} \\
        --output-dir ${prefix}_results \\
        ${contigs}

    tar -czf ${prefix}_antismash.tar.gz ${prefix}_results

    # To build the GFF3 file the scripts needs the regions.js file to be converted to json
    # In order to do that this process uses nodejs (using a patched version of the antismash container)

    echo ";var fs = require('fs'); fs.writeFileSync('./regions.json', JSON.stringify(recordData));" >> ${prefix}_results/regions.js

    node ${prefix}_results/regions.js

    mv regions.json ${prefix}_regions.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
        antiSMASH database: ${db_version}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_antismash.tar.gz
    touch ${prefix}_regions.json
    mkdir ${prefix}_results/
    touch ${prefix}_results/${prefix}.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
        antiSMASH database: ${db_version}
    END_VERSIONS
    """
}
