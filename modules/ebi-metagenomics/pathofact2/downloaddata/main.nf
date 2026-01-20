def downloadZenodoApiEntry(zenodo_id, max_retries = 3) {
    def api_url = "https://zenodo.org/api/records/${zenodo_id}"
    def conn = new URL(api_url).openConnection()
    conn.setRequestProperty('Accept', 'application/json')
    conn.setRequestProperty('User-Agent', "Nextflow ${nextflow.version ?: ''}".trim())
    conn.setConnectTimeout(10000)
    conn.setReadTimeout(30000)
    def api_text = conn.getInputStream().getText('UTF-8')
    def parser = new groovy.json.JsonSlurper()
    return parser.parseText(api_text)
}

process PATHOFACT2_DOWNLOADDATA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    output:
    path "pathofact_db/Models", emit: db
    tuple val("${task.process}"), val('pathofactdownloaddata'), eval("wget --version | head -n1 | cut -d' ' -f3"), topic: versions, emit: versions_pathofactdownloaddata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def zenodo_id = 18223764
    
    // Fetch API data to get version info
    def api_data = downloadZenodoApiEntry(zenodo_id)
    def db_version = api_data.metadata?.version ?: "unknown"
    
    // Use the direct download link from the files array
    // Zenodo provides multiple link options - we need the direct one
    def file_entry = api_data.files[0]
    def filename = file_entry.key
    def download_url = "https://zenodo.org/records/${zenodo_id}/files/${filename}"
    """
    # Using wget for reliable Zenodo downloads
    wget \
        ${args} \
        --no-verbose \
        --user-agent="Nextflow ${workflow.nextflow.version}" \
        --output-document=Models.tar.gz \
        "${download_url}"

    mkdir pathofact_db
    tar -xzf Models.tar.gz
    mv Models pathofact_db
    rm Models.tar.gz
    """

    stub:
    """
    mkdir -p pathofact_db/Models
    """
}
