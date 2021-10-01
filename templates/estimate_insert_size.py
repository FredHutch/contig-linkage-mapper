#!/usr/bin/env python3

import json
import gzip
import os
import pysam
import statistics
from collections import defaultdict

# Input data should be present in 'paired.bam'
fp = 'paired.bam'
assert os.path.exists(fp)

# Keep track of each read pair
read_pairs = defaultdict(dict)

# Open the input
with pysam.AlignmentFile(fp, "rb") as handle:

    # Iterate over each read
    for read in handle:

        # Assign the read pair label (R1/R2)
        read_pair = 'R1' if read.is_read1 else 'R2'

        # Make sure that we haven't seen an alignment for this read yet
        msg = "Found multiple alignments for %s from %s against the target"
        msg = msg % (read.query_name, "$specimen")
        assert read_pairs[read.query_name].get(read_pair) is None, msg

        # Assign the strand
        read_strand = 'rev' if read.is_reverse else 'fwd'

        # Get the position on the reference where the first position of the read aligned
        read_start = read.reference_end if read.is_reverse else read.reference_start

        # Save to the read_pairs dict
        read_pairs[read.query_name][read_pair] = dict(
            read_pair=read_pair,
            read_strand=read_strand,
            read_start=read_start
        )

# Now that we've read in all of the alignment information,
# go through and calculate the insert sizes

# Keep a running list of the insert sizes that we've found
insert_size_list = list()

# Iterate over the dict of dicts
for read_name, read_pair in read_pairs.items():

    # Both reads should have aligned
    msg = "Read was not aligned in a pair (%s / %s)" % (read_name, "$specimen")
    assert len(read_pair) == 2, msg

    # The reads should be on opposite strands
    msg = "Read pairs are aligned to the same strand"
    assert read_pair["R1"]["read_strand"] != read_pair["R2"]["read_strand"], msg

    # If R1 is aligned on the fwd strand
    if read_pair["R1"]["read_strand"] == "fwd":

        # Then the start position should be before R2
        msg = "Reads are not oriented inwards (%s / $specimen)" % read_name
        assert read_pair["R1"]["read_start"] < read_pair["R2"]["read_start"], msg

    # Otherwise
    else:

        # R1 is aligned on the reverse strand
        assert read_pair["R1"]["read_strand"] == "rev"

        # Then the start position should be after R2
        msg = "Reads are not oriented inwards (%s / $specimen)" % read_name
        assert read_pair["R1"]["read_start"] > read_pair["R2"]["read_start"], msg

    # Having passed all of the checks, add to the list of inserts
    insert_size_list.append(
        abs(read_pair["R1"]["read_start"] - read_pair["R2"]["read_start"])
    )
            
# Write out the statistics for the paired reads to a file
with open("insert_size.json", "w") as handle:

    json.dump(
        dict(
            n_pairs=len(insert_size_list),
            insert_median=statistics.median(insert_size_list),
            insert_mean=statistics.mean(insert_size_list),
            insert_stdev=statistics.stdev(insert_size_list)
        ),
        handle
    )