process MULTIQC {
    tag "$meta.id"
    label 'process_high_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    module 'multiqc/1.0dev'

    input:
    tuple val(meta), file ('fastqc/*')

    output:
    file "multiqc_report.html"
    file "multiqc_data"

    script:
    """
    multiqc .
    """
}