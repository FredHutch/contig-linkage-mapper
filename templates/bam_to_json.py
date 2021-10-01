#!/usr/bin/env python3

import json
import gzip
import os
import pysam

contigs_fp = "contigs.bam"
assert os.path.exists(contigs_fp)
target_fp = "target.bam"
assert os.path.exists(target_fp)

# Read a dict of reference lengths from the header
def read_reference_lengths(fp):
    assert os.path.exists(fp)
    with pysam.AlignmentFile(fp, "rb") as handle:
        return dict(zip(
            handle.references,
            handle.lengths
        ))

# Read in the position and orientation of each aligned read
def read_bam(fp):
    assert os.path.exists(fp)

    # Length dict
    length_dict = read_reference_lengths(fp)

    with pysam.AlignmentFile(fp, "rb") as handle:

        return [
            parse_read(read, length_dict)
            for read in handle
        ]

def parse_read(read, length_dict):
    return {
        "strand": 'rev' if read.is_reverse else 'fwd',
        "reference_name": read.reference_name,
        "reference_start": read.reference_end if read.is_reverse else read.reference_start,
        "reference_length": length_dict[read.reference_name],
        "query_name": read.query_name,
    }

def bam_to_json(bam_fp):

    json_fp = bam_fp.replace(".bam", ".json.gz")

    with gzip.open(json_fp, "wt") as handle:

        json.dump(read_bam(bam_fp), handle)

bam_to_json(target_fp)
bam_to_json(contigs_fp)
