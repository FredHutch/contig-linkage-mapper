#!/usr/bin/env python3

import numpy as np
import pandas as pd
import json

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

# Read in the summary of the read quality trimming
fastq_stats = json.load(open("cutadapt.json"))

n_reads = fastq_stats["read_counts"]["output"]
print("Total number of reads: %d" % n_reads)

# Read in the list of all junctions
junctions_df = pd.read_csv("junctions.csv.gz")
print(junctions_df)

def summarize_contig(d):
    # Count up the proportion of junctions which are 'fwd'
    prop_fwd = (d.orientation=='fwd').mean()

    return pd.Series(dict(
        pos_offset=d.pos_offset.median(),
        orientation='fwd' if prop_fwd >= 0.5 else 'rev',
        linked_contig_len=d.linked_contig_len.values[0],
        target_contig_len=d.target_contig_len.values[0],
        n_reads=d.shape[0],
        fpm=d.shape[0] / (n_reads / 1000000.)
    ))

# Summarize each linked contig
summary_df = junctions_df.groupby(
    "linked_contig",
).apply(
    summarize_contig
).reset_index()

summary_df.to_csv("summary.csv.gz", index=None)