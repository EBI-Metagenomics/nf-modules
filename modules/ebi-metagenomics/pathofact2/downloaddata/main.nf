process PATHOFACT2_DOWNLOADDATA {
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    output:
    path "pathofact_db/Models", emit: db
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget https://zenodo.org/records/14192463/files/DATABASES.tar.gz
    tar -xvf DATABASES.tar.gz
    mv DATABASES pathofact_db
    rm DATABASES.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS

    """
    stub:
    """
    mkdir -p pathofact_db/MODELS
    """
}
