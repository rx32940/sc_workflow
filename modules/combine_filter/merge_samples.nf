process MERGE_SAMPLES {
    tag "combine_filter"
    label 'process_high'
    publishDir "${params.outdir}/writes", mode: 'copy'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Rachel/envs/scanpy'

    input:
    val trigger  // Accepts the "true" value passed from ch_ready_to_merge

    output:
    path("merged_cts.h5ad")

    script:
    def samplesheet = file(params.samplesheet)
    def workdir = file(params.workdir ?: "${params.outdir}")

    def selected_batches_file = null

    // Case 1: user provides a comma-separated list
    if (params.selected_samples) {
        def selected_list = params.selected_samples.split(",").collect { it.trim() }.join("\n")
        selected_batches_file = file("${workDir}/selected_batches.txt")
        selected_batches_file.text = selected_list
    }

    // Case 2: user provides an actual file
    else if (params.selected_batches_file) {
        def candidate_file = file(params.selected_batches_file)
        if (candidate_file.exists()) {
            selected_batches_file = candidate_file
        }
    }

    // Prepare CLI argument only if needed
    def selected_flag = selected_batches_file ? "--selected_batches ${selected_batches_file}" : ""

    """

    # Run the main combine filter script
    /research/groups/northcgrp/home/common/Rachel/envs/scanpy/bin/python ${workflow.projectDir}/modules/combine_filter/merge_samples.py \
        --workdir ${workdir} \
        --samplesheet ${samplesheet} \
        --output merged_cts.h5ad \
        ${selected_flag}
    """
}
