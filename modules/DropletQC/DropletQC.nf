process DROPLET_QC {
    tag "$meta.id"
    label 'process_med'
    publishDir "${params.outdir}/DropletQC/${meta.id}", mode: 'copy'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Rachel/envs/R_nf_4.3.1'

    input:
    tuple val(meta), path(cranger_out)

    output:
    tuple val(meta), path("nf_ed_qc.csv"), path("empty_droplets_qc.png")

    script:
    """
    # Copy arc BAM files to expected names
    if [ -f ${cranger_out}/outs/gex_possorted_bam.bam ]; then
        cp ${cranger_out}/outs/gex_possorted_bam.bam ${cranger_out}/outs/possorted_genome_bam.bam
        cp ${cranger_out}/outs/gex_possorted_bam.bam.bai ${cranger_out}/outs/possorted_genome_bam.bam.bai
    fi

    Rscript ${workflow.projectDir}/modules/DropletQC/DropletQC_function.R \
    ${cranger_out}/outs \
    ${task.cpus}
    """
}