#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.target = false
params.contigs = false
params.reads = false
params.output = false
params.min_align_score = 50
params.min_qvalue = 20
params.cpus = 4
params.sigma = 1

// Set up the containers using params
params.container__pandas = "quay.io/fhcrc-microbiome/python-pandas:4a6179f"
params.container__bwa = "quay.io/fhcrc-microbiome/bwa:bwa.0.7.17__bcw.0.3.0I"
params.container__pysam = "quay.io/biocontainers/pysam:0.17.0--py36h61e5637_0"
params.container__cutadapt = "quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0"
params.container__multiqc = "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"

// Import the required processes
include { 
    quality_trim;
    bwa_index as bwa_index_target;
    bwa_index as bwa_index_contigs;
    parse_target;
    align_target;
    split_bam;
    estimate_insert_size;
    get_bam_headers;
    get_unmapped_reads;
    concat_contigs;
    align_contigs;
    bam_to_json;
    count_junctions;
    summarize_results;
    multiqc;
} from "./modules"

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/meta-integration-mapper <ARGUMENTS>

Required Arguments:
    --target        File containing target sequence in FASTA format
    --contigs       File(s) containing comparison contigs in FASTA format (gzip optional)
    --reads         File(s) containing paired-end WGS reads in FASTQ format (see below)
    --output        Path to output directory

Paired-end reads:
    WGS data must be specified in terms of file pairs, with syntax following the pattern
    described here: https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
    Note that this syntax will also determine the specimen names which are inferred.

Multiple input files:
    To specify more than one file for --contigs or --reads, wildcards can be used
    as well as listing multiple distinct paths separated by commmas.

Webpage: https://github.com/FredHutch/meta-integration-mapper
    """.stripIndent()
}


// Define the workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or does not provide the required inputs
    if (params.help || !params.target || !params.contigs || !params.reads || !params.output){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // Make a channel with all of the inputs reads
    reads_ch = Channel.fromFilePairs(params.reads)

    // Run cutadapt to trim low-quality bases from the input
    quality_trim(
        reads_ch
            .map {[it[0], it[1][0], it[1][1]]},
    )

    // Parse the target contig (and make sure there's only one contig in that file)
    parse_target(file(params.target))

    // Index the target
    bwa_index_target(parse_target.out)

    // Align all input reads to the target contig
    align_target(
        // Join each pair of FASTQ files with the contig FASTA
        quality_trim.out.bam,
        bwa_index_target.out
    )

    // Split up the aligned reads into paired and unpaired sequences
    split_bam(
        align_target.out.bam
    )

    // Estimate the insert size from the paired alignments to the target contig
    estimate_insert_size(
        split_bam.out.paired
    )

    // Extract the headers of R1 and R2 which were aligned as single-ended
    get_bam_headers(
        split_bam.out.unpaired
    )

    // Get the other half of the single-end aligned reads
    get_unmapped_reads(
        get_bam_headers.out
            .join(quality_trim.out.bam)
    )

    // Merge all of the contigs
    concat_contigs(
        Channel.fromPath(params.contigs.tokenize(",")).toSortedList()
    )

    // Index the merged contigs
    bwa_index_contigs(
        concat_contigs.out
    )

    // Align the subset of reads which match only one end to the target contig
    align_contigs(
        get_unmapped_reads.out,
        bwa_index_contigs.out
    )

    // Convert the alignment information to JSON for easier parsing downstream
    bam_to_json(
        align_contigs.out.bam
            .join(align_target.out.bam)
    )

    // Count up the junction points
    count_junctions(
        bam_to_json.out
            .join(estimate_insert_size.out)
    )

    // Summarize all of the results
    summarize_results(
        count_junctions.out
            .join(quality_trim.out.log)
    )

    // Generate an overall QC report
    multiqc(
        align_target.out.log.toSortedList(),       
        align_contigs.out.log.toSortedList(),       
    )

}