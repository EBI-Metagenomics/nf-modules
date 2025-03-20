
process CONVERSTCMSCANTOCMSEARCH {
    tag "$meta.id"
    label 'process_single'

    container "microbiome-informatics/convert_cmscan_to_cmsearch:1"

    input:
    tuple val(meta), path(cmscan_tblout)

    output:
    tuple val(meta), path("*cmsearch.tbl"), emit: cmsearch_tblout

    script:

    """
    convert-cmscan-to-cmsearch-tblout.py --input $cmscan_tblout --output ${meta.id}.cmsearch.tbl
    """
}
