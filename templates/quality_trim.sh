#!/bin/bash

set -Eeuo pipefail

cutadapt \
    --pair-filter=any \
    --quality-cutoff=${params.min_qvalue} \
    --minimum-length=${params.min_align_score} \
    -o trimmed.R1.fastq.gz \
    -p trimmed.R1.fastq.gz \
    --json="${specimen}.cutadapt.json" \
    $R1 \
    $R2
