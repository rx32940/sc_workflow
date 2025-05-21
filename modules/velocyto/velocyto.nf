process VELOCYTO_RUN10X {
    tag "$meta.id"
    label 'process_med'
    publishDir "${params.outdir}/Velocyto", mode: 'copy'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Rachel/envs/velocyto'

    input:
    tuple val(meta), path(cranger_out)

    output:
    tuple val(meta), path("${meta.id}.loom")

    script:
    def gtf_file = meta.gtf

    """
    # Sort the BAM file by cell barcode
    samtools sort -t CB -O BAM -@ ${task.cpus} \
        -o ${cranger_out}/outs/cellsorted_possorted_genome_bam.bam \
        ${cranger_out}/outs/possorted_genome_bam.bam

    # Run Velocyto
    velocyto run10x ${cranger_out} ${gtf_file}
    
    # Move the loom file to the final directory
    mv ${cranger_out}/velocyto/*.loom ${meta.id}.loom
    """
}
