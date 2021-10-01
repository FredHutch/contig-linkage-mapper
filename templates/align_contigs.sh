#!/bin/bash

set -Eeuo pipefail

echo "Aligning WGS data from $specimen"
ls -lahtr

# Align and keep all alignments (-a),
# remove all unmapped reads, and
# write to ${specimen}.bam
bwa mem \
    -t ${task.cpus} \
    -T ${params.min_align_score} \
    -a \
    ref.fasta \
    ${fastq} \
    2>${specimen}.contigs.log | \
    samtools \
        view \
        -h \
        -b \
        -F 4 \
    > aligned.bam
