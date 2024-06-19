process INTERPROSCAN {
    tag "$meta.id"
    label 'process_long'

    container 'microbiome-informatics/genomes-pipeline.ips:5.62-94.0'
    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return "--bind ${interproscan_db}/data:/opt/interproscan-5.62-94.0/data"
        } else {
            return "-v ${interproscan_db}/data:/opt/interproscan-5.62-94.0/data"
        }
    }

    input:
    tuple val(meta), path(fasta)
    tuple path(interproscan_db), val(db_version)
    val(out_ext)

    output:
    tuple val(meta), path('*.tsv') , optional: true, emit: tsv
    tuple val(meta), path('*.xml') , optional: true, emit: xml
    tuple val(meta), path('*.gff3'), optional: true, emit: gff3
    tuple val(meta), path('*.json'), optional: true, emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def appl = "-appl TIGRFAM,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,PfamA,MobiDBLite"
    if ( args.contains("-appl") ) {
        appl = ""
    }
    switch ( out_ext ) {
        case "tsv": break
        case "xml": break
        case "gff3": break
        case "json": break
        default:
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting tsv as fallback");
            break
    }

    //  -dp (disable precalculation) is on so no online dependency
    """
    interproscan.sh \\
        -cpu $task.cpus \\
        -i $fasta \\
        -f ${out_ext} \\
        -dp \\
        ${appl} \\
        ${args} \\
        -o ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    switch ( out_ext ) {
        case "tsv": break
        case "xml": break
        case "gff3": break
        case "json": break
        default:
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting tsv as fallback");
            break
    }

    """
    touch ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """
}
