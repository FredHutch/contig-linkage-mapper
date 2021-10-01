#!/usr/bin/env python3

import pysam
import gzip
import os

specimen = "$specimen"
print("Processing %s" % specimen)

R1_txt = "R1.txt.gz"
assert os.path.exists(R1_txt)
R2_txt = "R2.txt.gz"
assert os.path.exists(R2_txt)
R1_fastq = "R1.fastq.gz"
assert os.path.exists(R1_fastq)
R2_fastq = "R2.fastq.gz"
assert os.path.exists(R2_fastq)

# Function to read in a gzipped text file
def read_txt_gz(fp):
    return set([line.rstrip("\\n") for line in gzip.open(fp, "rt")])

# Read in the names of all R1's whose mates were unmapped
R1_mates_unmapped = read_txt_gz(R1_txt)
print("R1's with unmapped mates: %d" % len(R1_mates_unmapped))

# Read in the names of all R2's whose mates were unmapped
R2_mates_unmapped = read_txt_gz(R2_txt)
print("R2's with unmapped mates: %d" % len(R2_mates_unmapped))

# Function to extract sequences with a particular header
def get_fastq_headers(fp, headers):

    # Open the file
    with pysam.FastxFile(fp) as handle:
    
        # Iterate over the sequences and return all matches
        return [
            entry
            for entry in handle
            if entry.name in headers
        ]

# Get the R2's which were unmapped, but whose mates were mapped
R1_mates = get_fastq_headers(R2_fastq, R1_mates_unmapped)
print("Unmapped R2's with mapped mates: %d" % len(R1_mates))

# Get the R1's which were unmapped, but whose mates were mapped
R2_mates = get_fastq_headers(R1_fastq, R2_mates_unmapped)
print("Unmapped R1's with mapped mates: %d" % len(R2_mates))

# Make sure that the numbers add up
assert len(R1_mates) == len(R1_mates_unmapped)
assert len(R2_mates) == len(R2_mates_unmapped)

# Write out to a file
with gzip.open("unmapped.fastq.gz", "wt") as handle:
    for entry in R1_mates + R2_mates:
        handle.write(str(entry) + "\\n")
