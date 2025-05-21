process CELLBENDER_REMOVE_BACKGROUND {
    tag "$meta.id"
    label 'process_gpu_high'
    publishDir "${params.outdir}/cellbender", mode: 'copy'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Rachel/envs/cellbender'

    input:
    tuple val(meta), path(cranger_out)

    output:
    tuple val(meta), path("${meta.id}_cb")

    """
    mkdir -p ${meta.id}_cb
    
    cellbender remove-background \
        --cuda \
        --input ${cranger_out}/outs/raw_feature_bc_matrix.h5 \
        --output ${meta.id}_cb/cellbender_output.h5
    """
}
