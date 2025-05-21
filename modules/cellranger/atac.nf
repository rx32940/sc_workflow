process CELLRANGER_ATAC_COUNT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/cellranger-atac", mode: 'copy'

    module 'cellranger-atac/2.1.0'

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("${meta.id}")

    script:
    def prefix = meta.id
    def genome_ref = meta.genome_ref

    """

    cellranger-atac count \
      --id ${prefix} \
      --sample ${prefix} \
      --fastqs ${fastq_dir} \
      --reference ${genome_ref} \
      --localcores ${task.cpus} \
      --localmem ${task.memory.toGiga()} 
    """
}