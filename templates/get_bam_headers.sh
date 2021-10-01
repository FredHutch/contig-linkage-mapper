#!/bin/bash

# Get the headers of all R1's
samtools view $bam -f 64 | cut -f 1 | gzip -c > "R1.txt.gz"

# Get the headers of all R2's
samtools view $bam -f 128 | cut -f 1 | gzip -c > "R2.txt.gz"
