#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Make sure that the target file only contains a single sequence
process parse_target {
    container "${params.container__pandas}"
    cpus 1
    
    input:
    path contig_fasta

    output:
    path contig_fasta

    script:
    template 'parse_target.py'

}

// Perform quality trimming on the input FASTQ data
process quality_trim {
    container "${params.container__cutadapt}"
    cpus 1
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path(R1), path(R2), emit: bam
    tuple val(specimen), path("${specimen}.cutadapt.json"), emit: log

    script:
    template 'quality_trim.sh'

}

// Index a FASTA file for alignment with BWA
process bwa_index {
    container "${params.container__bwa}"
    cpus 1
    
    input:
    path "ref.fasta"

    output:
    path "ref.fasta*"

    script:
    """
    bwa index ref.fasta
    """

}

// Align paired-end reads against the target contig (paired-end)
process align_target {
    container "${params.container__bwa}"
    cpus params.cpus
    
    input:
    tuple val(specimen), path(R1), path(R2)
    path "*"

    output:
    tuple val(specimen), path("aligned.bam"), emit: bam
    path "${specimen}.target.log", emit: log

    script:
    template 'align_target.sh'

}

// Align paired-end reads against the other contigs (single-end)
process align_contigs {
    container "${params.container__bwa}"
    cpus params.cpus
    
    input:
    tuple val(specimen), path(fastq)
    path "*"

    output:
    tuple val(specimen), path("aligned.bam"), emit: bam
    path "${specimen}.contigs.log", emit: log

    script:
    template 'align_contigs.sh'

}

// Split up a BAM into reads with are (a) both mapped and (b) only 1/2 mapped
process split_bam {
    container "${params.container__bwa}"
    cpus 1
    
    input:
    tuple val(specimen), path(bam)

    output:
    tuple val(specimen), path("paired.bam"), emit: paired
    tuple val(specimen), path("unpaired.bam"), emit: unpaired

    script:
    template 'split_bam.sh'

}

// Estimate the insert size for a specimen based on paired alignments to the target
process estimate_insert_size {
    container "${params.container__pysam}"
    cpus 1
    
    input:
    tuple val(specimen), path("paired.bam")

    output:
    tuple val(specimen), path("insert_size.json")

    script:
    template 'estimate_insert_size.py'

}

// Get the FASTQ headers for both R1 and R2 from a BAM
process get_bam_headers {
    container "${params.container__bwa}"
    cpus 1
    
    input:
    tuple val(specimen), path(bam)

    output:
    tuple val(specimen), path("R1.txt.gz"), path("R2.txt.gz")

    script:
    template 'get_bam_headers.sh'

}

// Get the FASTQ sequences of reads which were unmapped
process get_unmapped_reads {
    container "${params.container__pysam}"
    cpus 1
    
    input:
    tuple val(specimen), path("R1.txt.gz"), path("R2.txt.gz"), path("R1.fastq.gz"), path("R2.fastq.gz")

    output:
    tuple val(specimen), path("unmapped.fastq.gz")

    script:
    template 'get_unmapped_reads.py'

}

// Combine a set of FASTA sequences
process concat_contigs {
    container "${params.container__bwa}"
    cpus 1
    
    input:
    path "input/*"

    output:
    path "contigs.fasta.gz"

    """#!/bin/bash

    gunzip -c input/* | gzip -c > contigs.fasta.gz
    """

}

// Convert the alignment information to JSON
process bam_to_json {
    container "${params.container__pysam}"
    cpus 1
    
    input:
    tuple val(specimen), path("contigs.bam"), path("target.bam")

    output:
    tuple val(specimen), path("contigs.json.gz"), path("target.json.gz")
    
    script:
    template 'bam_to_json.py'

}

// Count up junction points
process count_junctions {
    container "${params.container__pandas}"
    cpus 1
    
    input:
    tuple val(specimen), path("contigs.json.gz"), path("target.json.gz"), path("insert_size.json")

    output:
    tuple val(specimen), path("junctions.csv.gz"), optional: true
    
    script:
    template 'count_junctions.py'

}

// Summarize the number of junctions per contig, per specimen
process summarize_results {
    container "${params.container__pandas}"
    publishDir "${params.output}", mode: 'copy', overwrite: true
    cpus 1
    
    input:
    tuple val(specimen), path("junctions.csv.gz"), path("cutadapt.json")

    output:
    path "${specimen}.csv.gz"
    
    script:
    template 'summarize_results.py'

}

// Generate a quality report for each specimen
process multiqc {
    container "${params.container__multiqc}"
    cpus 1
    
    input:
    path "target_alignments/*"
    path "contig_alignments/*"

    output:
    path "multiqc.html"
    
    script:
    template 'multiqc.sh'

}
