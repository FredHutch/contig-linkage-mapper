#!/usr/bin/env python3

import gzip
import os

# Get the input file path
fp = "${contig_fasta}"

# Make sure that the file exists
assert os.path.exists(fp), "Input file not found"

# If the file ends with ".gz"
if fp.endswith(".gz"):

    # Open it with gzip
    handle = gzip.open(fp, "rt")

# If it does not end with ".gz"
else:

    # Open it as a regular file
    handle = open(fp, "r")

# Count up the number of lines
nlines = 0
# Count up the number of headers
nheaders = 0

# Iterate over each line
for line in handle:

    # If the line is a header
    if line.startswith(">"):
        nheaders += 1

    # Count the line
    nlines += 1

assert nlines > 0, "Target file must contain a single contig"
assert nheaders > 0, "Target file must contain a FASTA header starting with '>'"
assert nheaders == 1, "Target file may only contain a single sequence"
assert nlines > 1, "Target file must contain both a header and nucleotide sequence"

print("Number of lines: %s" % nlines)
print("Number of headers: %s" % nheaders)

# Close the file handle
handle.close()
