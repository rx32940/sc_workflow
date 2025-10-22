#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ }                    from './modules/cat/fastq'
include { CELLRANGER_COUNT }             from './modules/cellranger/count'
include { CELLRANGER_ATAC_COUNT }        from './modules/cellranger/atac'
include { CELLBENDER_REMOVE_BACKGROUND } from './modules/cellbender/cellbender'
include { VELOCYTO_RUN10X }              from './modules/velocyto/velocyto'
include { DROPLET_QC }                   from './modules/DropletQC/DropletQC'
include { MERGE_SAMPLES }                from './modules/combine_filter/merge_samples'

workflow {

    if(params.modality == "atac") {
    
        // Set the genome reference based on species or user-supplied path
        def human_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Human/Ensembl/CellRanger/GRCh38/refdata-cellranger-arc-GRCh38-2024-A"
        def mouse_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Mouse/Ensembl/CellRanger/GRCm39/refdata-cellranger-arc-GRCm39-2024-A"

        def genome_ref = params.genome_ref ? file(params.genome_ref).toString() :
                        (params.species == "mouse" ? mouse_ref : human_ref)
        
        if (!file(genome_ref).exists()) {
            error "Genome reference not found: ${genome_ref}"
        }
        
        // Set the GTF path
        def gtf_file = params.gtf ? file(params.gtf).toString() :
                        file("${genome_ref}/genes/genes.gtf").toString()
        
        if (!file(gtf_file).exists()) {
            error "GTF file not found at expected location: ${gtf_file}"
        }
        
        // Read the samplesheet
        Channel
            .fromPath(params.samplesheet)
            .splitCsv(skip: 1, sep: ',', strip: true)
            .map { id, fastq_dir_str ->
                def fastq_dir = file(fastq_dir_str)
        
                // Create metadata
                def meta = [id: id, genome_ref: genome_ref, gtf: gtf_file]
        
                // Check FASTQ presence
                def fastqs = fastq_dir.listFiles()?.findAll { it.name.endsWith(".fastq.gz") } ?: []
                if (fastqs.isEmpty()) {
                    error "No FASTQ files found in directory: ${fastq_dir}"
                }
        
                return [meta, fastq_dir]
            }
            .set { ch_fastq }


      CELLRANGER_ATAC_COUNT(ch_fastq)

    } else {

    // Set the genome reference based on species or user-supplied path
    def human_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Human/Ensembl/CellRanger/GRCh38/refdata-gex-GRCh38-2024-A"
    def mouse_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Mouse/Ensembl/CellRanger/GRCm39/refdata-gex-GRCm39-2024-A"

    def genome_ref = params.genome_ref ? file(params.genome_ref).toString() :
                       (params.species == "mouse" ? mouse_ref : human_ref)

    if (!file(genome_ref).exists()) {
        error "Genome reference not found: ${genome_ref}"
    }

    // Set the GTF path
    def gtf_file = params.gtf ? file(params.gtf).toString() :
                    file("${genome_ref}/genes/genes.gtf").toString()

    if (!file(gtf_file).exists()) {
        error "GTF file not found at expected location: ${gtf_file}"
    }

    // Read the samplesheet
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(skip: 1, sep: ',', strip: true)
        .map { id, fastq_dir_str ->
            def fastq_dir = file(fastq_dir_str)

            // Create metadata
            def meta = [id: id, genome_ref: genome_ref, gtf: gtf_file]

            // Check FASTQ presence
            def fastqs = fastq_dir.listFiles()?.findAll { it.name.endsWith(".fastq.gz") } ?: []
            if (fastqs.isEmpty()) {
                error "No FASTQ files found in directory: ${fastq_dir}"
            }

            return [meta, fastq_dir]
        }
        .set { ch_fastq }

    // Run Cell Ranger
    CELLRANGER_COUNT(ch_fastq)
        .set { ch_cellranger_out }

    // Trigger post-count processes
    ch_cellranger_out
        .map { meta, outs_dir -> [meta, outs_dir] }
        .set { ch_post_count }

    CELLBENDER_REMOVE_BACKGROUND(ch_post_count)
        .set { ch_cellbender_done }

    VELOCYTO_RUN10X(ch_post_count)
        .set { ch_velocyto_done }

    DROPLET_QC(ch_post_count)
        .map { meta, csv, _ -> tuple(meta, csv) }  // ignore the .png
        .set { ch_dropletcq_done }


    // Sync signal after all three finish
    ch_cellbender_done
        .combine(ch_velocyto_done)
        .combine(ch_dropletcq_done)
        .collect()
        .map { _ -> true }
        .set { ch_ready_to_merge }

    // Run MERGE_SAMPLES when all is done
    MERGE_SAMPLES(ch_ready_to_merge)
    }
}
