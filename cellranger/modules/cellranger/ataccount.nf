// https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count
process CELLRANGER_ATAC_COUNT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/scatac", mode: 'copy'

    module 'cellranger-atac/2.1.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("**/outs/**")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir fastq_all
    mv *.fastq.gz fastq_all
    cellranger-atac count \
      --id ${prefix} \
      --sample ${prefix} \
      --fastqs fastq_all \
      --reference ${params.atac_ref} \
      --localcores ${task.cpus} \
      --localmem ${task.memory.toGiga()} 
    """
}