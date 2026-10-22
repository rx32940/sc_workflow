process RHAPSODY_TO_CELLRANGER {
    tag "$meta.id"
    label 'process_high'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Rachel/envs/velocyto'

    input:
    tuple val(meta), path(bam), path(unfiltered)

    output:
    tuple val(meta), path("${meta.id}")

    script:
    // meta.filtered_h5 (optional) is a user-supplied filtered h5; otherwise filter by params.min_umi
    def filtered_arg = meta.filtered_h5 ? "--filtered ${meta.filtered_h5}" : "--min_umi ${params.min_umi}"
    """
    python ${workflow.projectDir}/modules/rhapsody/rhapsody_to_cellranger.py \
        --bam ${bam} \
        --unfiltered ${unfiltered} \
        ${filtered_arg} \
        --outdir ${meta.id} \
        --genome ${params.species == 'mouse' ? 'GRCm39' : 'GRCh38'} \
        --threads ${task.cpus}
    """
}
