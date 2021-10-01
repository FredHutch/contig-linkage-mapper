#!/bin/bash

set -Eeuo pipefail

echo """
Splitting up reads in $bam

"""
ls -lahtr

# Reads which have mapped mates
samtools \
    view \
    -h \
    -b \
    -F 8 \
    $bam > "paired.bam"

# Reads which have UNmapped mates
samtools \
    view \
    -h \
    -b \
    -f 8 \
    $bam > "unpaired.bam"

echo """
DONE

"""

ls -lahtr
