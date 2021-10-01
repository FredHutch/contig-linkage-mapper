#!/bin/bash

set -Eeuo pipefail

echo "Aligning WGS data from $specimen"
ls -lahtr

# Align, remove all unmapped reads, and write to ${specimen}.bam
bwa mem \
    -t ${task.cpus} \
    -T ${params.min_align_score} \
    ref.fasta \
    ${R1} \
    ${R2} \
    2>${specimen}.target.log | \
    samtools \
        view \
        -h \
        -b \
        -F 4 \
    > aligned.bam
