process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/cellranger", mode: 'copy'

    module 'cellranger/8.0.1'

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("${meta.id}")

    script:
    def prefix = meta.id
    def genome_ref = meta.genome_ref

    """
    # Run cellranger count with optimized resource management
    cellranger count \
      --id ${prefix} \
      --sample ${prefix} \
      --fastqs ${fastq_dir} \
      --transcriptome ${genome_ref} \
      --localcores ${task.cpus} \
      --localmem ${task.memory.toGiga()} \
      --create-bam true \
      --disable-ui \
      --chemistry ${params.chemistry ?: 'auto'}
    """
}
