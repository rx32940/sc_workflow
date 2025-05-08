#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ } from './modules/cat/fastq'
include { FASTQC } from './modules/fastqc/fastqc'
include { MULTIQC } from './modules/fastqc/multiqc'
include { CELLRANGER_RNA_COUNT } from './modules/cellranger/rnacount'
include { CELLRANGER_ATAC_COUNT } from './modules/cellranger/ataccount'

workflow {
  Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map {
      row ->
        def fastq = row.findAll { it.key != 'sample' }
        return [[id:row.sample], [fastq.values()]]
    }
    .groupTuple(by: [0])
    .map {meta, fastq -> [meta, fastq.flatten()] }
    .set { ch_fastq }

    if(params.modality == "atac") {
      CELLRANGER_ATAC_COUNT(ch_fastq)
    } else {
      CELLRANGER_RNA_COUNT(ch_fastq)
    }

    // not yet used
    //ch_cat_fastq | FASTQC | collect |  MULTIQC

}
