#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ }                    from './modules/cat/fastq'
include { CELLRANGER_COUNT }             from './modules/cellranger/count'
include { CELLRANGER_ARC_COUNT }        from './modules/cellranger/arc'
include { CELLRANGER_ATAC_COUNT }        from './modules/cellranger/atac'
include { CELLBENDER_REMOVE_BACKGROUND } from './modules/cellbender/cellbender'
include { VELOCYTO_RUN10X }              from './modules/velocyto/velocyto'
include { DROPLET_QC }                   from './modules/DropletQC/DropletQC'
include { MERGE_SAMPLES }                from './modules/combine_filter/merge_samples'
include { RHAPSODY_TO_CELLRANGER }       from './modules/rhapsody/rhapsody_to_cellranger'

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
        // ARC, COUNT or RHAPSODY: run mapping (or convert existing BD Rhapsody output) + downstream analysis
        
        if(params.modality == "arc") {
            
            def human_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Human/Ensembl/CellRanger/GRCh38/refdata-cellranger-arc-GRCh38-2024-A"
            def mouse_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Mouse/Ensembl/CellRanger/GRCm39/refdata-cellranger-arc-GRCm39-2024-A"

            def genome_ref = params.genome_ref ? file(params.genome_ref).toString() :
                            (params.species == "mouse" ? mouse_ref : human_ref)

            if (!file(genome_ref).exists()) {
                error "Genome reference not found: ${genome_ref}"
            }

            Channel
                .fromPath(params.samplesheet)
                .splitCsv(header: true, sep: ',', strip: true)
                .map { row -> 
                    tuple(row.sample, row.fastqs, row.library_type)
                }
                .groupTuple()
                .map { sample, fastq_dirs, library_types ->
                    def gtf_file = params.gtf ? file(params.gtf).toString() :
                                    file("${genome_ref}/genes/genes.gtf").toString()
                    def meta = [id: sample, genome_ref: genome_ref, gtf: gtf_file, modality: "arc"]
                    def libraries_content = "fastqs,sample,library_type\n"
                    [fastq_dirs, library_types].transpose().each { dir, type ->
                        libraries_content += "${dir},${sample},${type}\n"
                    }
                    return [meta, libraries_content, fastq_dirs]
                }
                .set { ch_libraries }

            CELLRANGER_ARC_COUNT(ch_libraries)
                .set { ch_cellranger_out }
                

        } else if(params.modality == "rhapsody") {
            // BD Rhapsody output already mapped by the Rhapsody pipeline: convert it to a
            // Cell Ranger-like layout instead of running Cell Ranger
            // samplesheet: id,rhapsody_dir[,raw_h5,filtered_h5]
            //   raw_h5 / filtered_h5 are optional 10x-style h5 files you made yourself (barcodes must be the
            //   Rhapsody cell labels in the BAM). Without raw_h5 the Rhapsody *_Unfiltered_MEX.zip is converted;
            //   without filtered_h5 the raw matrix is filtered to barcodes with >= params.min_umi UMIs.

            def human_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Human/Ensembl/CellRanger/GRCh38/refdata-gex-GRCh38-2024-A"
            def mouse_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Mouse/Ensembl/CellRanger/GRCm39/refdata-gex-GRCm39-2024-A"

            // GTF used by velocyto; defaults to the Cell Ranger GENCODE GTF, override with --gtf
            // to use the annotation matching your Rhapsody reference
            def genome_ref = params.genome_ref ? file(params.genome_ref).toString() :
                           (params.species == "mouse" ? mouse_ref : human_ref)
            def gtf_file = params.gtf ? file(params.gtf).toString() :
                            file("${genome_ref}/genes/genes.gtf").toString()

            if (!file(gtf_file).exists()) {
                error "GTF file not found: ${gtf_file}"
            }

            Channel
                .fromPath(params.samplesheet)
                .splitCsv(header: true, sep: ',', strip: true)
                .map { row ->
                    def rdir = file(row.rhapsody_dir)
                    def pick = { String pattern ->
                        def hits = rdir.listFiles()?.findAll { it.name ==~ pattern } ?: []
                        if (hits.size() != 1) {
                            error "Expected exactly one file matching '${pattern}' in ${rdir}, found ${hits.size()}"
                        }
                        return hits[0]
                    }
                    def user_h5 = { String v ->
                        if (!file(v).exists()) { error "h5 file not found: ${v}" }
                        return file(v)
                    }
                    def meta = [id: row.id, genome_ref: genome_ref, gtf: gtf_file, modality: "count",
                                filtered_h5: row.filtered_h5 ? user_h5(row.filtered_h5).toString() : null]
                    return [meta,
                            pick(/.*_Bioproduct\.bam/),
                            row.raw_h5 ? user_h5(row.raw_h5) : pick(/.*_RSEC_MolsPerCell_Unfiltered_MEX\.zip/)]
                }
                .set { ch_rhapsody }

            RHAPSODY_TO_CELLRANGER(ch_rhapsody)
                .set { ch_cellranger_out }

        } else {
            // COUNT modality
            
            def human_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Human/Ensembl/CellRanger/GRCh38/refdata-gex-GRCh38-2024-A"
            def mouse_ref = "/research/groups/northcgrp/projects/Northcott_Bioinformatics/northcgrp/References/Mouse/Ensembl/CellRanger/GRCm39/refdata-gex-GRCm39-2024-A"

            def genome_ref = params.genome_ref ? file(params.genome_ref).toString() :
                           (params.species == "mouse" ? mouse_ref : human_ref)

            if (!file(genome_ref).exists()) {
                error "Genome reference not found: ${genome_ref}"
            }

            def gtf_file = params.gtf ? file(params.gtf).toString() :
                            file("${genome_ref}/genes/genes.gtf").toString()

            if (!file(gtf_file).exists()) {
                error "GTF file not found at expected location: ${gtf_file}"
            }

            Channel
                .fromPath(params.samplesheet)
                .splitCsv(skip: 1, sep: ',', strip: true)
                .map { id, fastq_dir_str ->
                    def fastq_dir = file(fastq_dir_str)
                    def meta = [id: id, genome_ref: genome_ref, gtf: gtf_file, modality: "count"]
                    def fastqs = fastq_dir.listFiles()?.findAll { it.name.endsWith(".fastq.gz") } ?: []
                    if (fastqs.isEmpty()) {
                        error "No FASTQ files found in directory: ${fastq_dir}"
                    }
                    return [meta, fastq_dir]
                }
                .set { ch_fastq }

            CELLRANGER_COUNT(ch_fastq)
                .set { ch_cellranger_out }
        }

        // Downstream analysis (shared for both arc and count)
        ch_cellranger_out
            .map { meta, outs_dir -> [meta, outs_dir] }
            .set { ch_post_count }

        CELLBENDER_REMOVE_BACKGROUND(ch_post_count)
            .set { ch_cellbender_done }

        VELOCYTO_RUN10X(ch_post_count)
            .set { ch_velocyto_done }

        DROPLET_QC(ch_post_count)
            .map { meta, csv, _ -> tuple(meta, csv) }
            .set { ch_dropletcq_done }

        ch_cellbender_done
            .combine(ch_velocyto_done)
            .combine(ch_dropletcq_done)
            .collect()
            .map { _ -> true }
            .set { ch_ready_to_merge }

        MERGE_SAMPLES(ch_ready_to_merge)
    }
}