process CELLRANGER_ARC_COUNT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/cellranger", mode: 'copy'

    module 'cellranger-arc/2.0.2'

    input:
    tuple val(meta), val(libraries_csv), val(fastq_dirs)

    output:
    tuple val(meta), path("${meta.id}")

    script:
    def prefix = meta.id
    def genome_ref = meta.genome_ref

    """
    # Create libraries.csv
cat > libraries.csv << 'EOF'
${libraries_csv}
EOF

    cellranger-arc count \
    --id ${prefix} \
    --reference ${genome_ref} \
    --libraries libraries.csv \
    --localcores ${task.cpus} \
    --localmem ${task.memory.toGiga()}
    """
}
