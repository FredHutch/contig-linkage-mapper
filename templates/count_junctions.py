#!/usr/bin/env python3

import gzip
import json
import os
import pandas as pd

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

# Get the sigma parameter
sigma = float("$params.sigma")
print("Using sigma=%s" % sigma)

# Read in the summary of estimated insert size for this specimen
insert_size_fp = "insert_size.json"
assert os.path.exists(insert_size_fp)
with open(insert_size_fp, "r") as handle:
    insert_size = json.load(handle)

# Calculate the distance threshold for 'end' alignments
end_threshold = int(insert_size["insert_mean"] + (insert_size["insert_stdev"] * sigma))
print("Using end threshold of %sbp" % end_threshold)

contigs_fp = "contigs.json.gz"
assert os.path.exists(contigs_fp)
target_fp = "target.json.gz"
assert os.path.exists(target_fp)

def read_json(fp):
    assert os.path.exists(fp)
    with gzip.open(fp, "rt") as handle:
        return pd.DataFrame(
            json.load(handle)
        )

# For each set of alignments, calculate the distance to
# the end of the contig that it is pointing towards, and
# whether that classifies it as 'internal' or 'end' based
# on the --sigma threshold provided by the user
def classify_junction_position(df):

    # Calculate the distance to the end of the contig that
    # the read is pointing towards
    df = df.assign(
        distance_to_end=df.apply(
            calc_distance_to_end,
            axis=1
        )
    )

    # Classify each read as 'internal' or 'end' based on the
    # distance to the end of the contig
    df = df.assign(
        position=df.apply(
            lambda r: 'end' if r.distance_to_end <= end_threshold else 'internal',
            axis=1
        )
    )

    return df


def calc_distance_to_end(r):

    # If the target read is aligned in the forward direction
    if r['strand'] == 'fwd':

        # Then we use the distance to the 3' end of the contig to determine
        # whether the position is 'internal' or 'end'
        return r['reference_length'] - r['reference_start']

    # Otherwise
    else:

        # The orientation must be in the reverse direction
        assert r['strand'] == 'rev'

        # Then we use the distance to the 5' end of the contig to determine
        # whether the position is 'internal' or 'end'
        return r['reference_start']


# Function to calculate the positional offset of the linked contig
def calc_positional_offset(target_aln, contig_aln):

    # The components which may be needed are
    target_start = target_aln['reference_start']
    insert_length = int(insert_size['insert_mean'])
    linked_start = contig_aln['reference_start']
    linked_length = contig_aln['reference_length']

    # If the two contigs are oriented in the same direction
    if target_aln['strand'] != contig_aln['strand']:

        # If the read on the target contig is oriented in the forward direction
        if target_aln['strand'] == 'fwd':

            # The positional offset is the position on the target
            # contig, plus the insert length, minus the position
            # on the linked contig
            return target_start + insert_length - linked_start

        # If the read on the target contig is oriented in the reverse direction
        else:

            assert target_aln['strand'] == 'rev'

            # The positional offset is the position on the target
            # contig, minus the insert length, minus the position
            # on the linked contig
            return target_start - insert_length - linked_start

    # If the two contigs are oriented in the same direction
    # (relative to the contig that they have been aligned to)
    else:

        # If the read on the target contig is oriented in the forward direction
        if target_aln['strand'] == 'fwd':

            # The positional offset is the position on the target
            # contig, plus the insert length, plus the position
            # on the linked contig, minus the length of the linked contig
            return target_start + insert_length + linked_start - linked_length

        # If the read on the target contig is oriented in the reverse direction
        else:

            assert target_aln['strand'] == 'rev'

            # The positional offset is the position on the target
            # contig, minus the insert length, plus the position
            # on the linked contig, minus the length of the linked contig
            return target_start - insert_length + linked_start - linked_length

# Read in the data and process it
contigs_df = classify_junction_position(read_json(contigs_fp))
target_df = classify_junction_position(read_json(target_fp))

# Keep a list of junctions
junction_list = list()

# Iterate over each read which aligned to the target
for read_name, target_aln in target_df.set_index("query_name").iterrows():

    # Find the pair for this read in the contig alignment
    ix = (contigs_df['query_name'] == read_name)

    # If the other half of this read did not align
    if ix.sum() == 0:
        
        # Skip it
        continue

    # Iterate over each alignment position across the contigs
    for _, contig_aln in contigs_df.loc[ix].iterrows():

        # Figure out the relative orientation of the linked contig
        orientation = 'fwd' if target_aln['strand'] == contig_aln['strand'] else 'rev'

        # Figure out the positional offset of the linked contig
        pos_offset = calc_positional_offset(target_aln, contig_aln)

        # Keep track of the data for this junction
        junction_list.append(
            dict(
                # The direction of the linkage is given by the strand of the
                # read aligning to the target contig
                direction=target_aln['strand'],
                # The relative orientation of the linked contig
                orientation=orientation,
                # Position of the read in the target contig
                target_position=target_aln['position'],
                # First aligned base of the fragment on the target contig
                target_start=target_aln['reference_start'],
                # Position of the read in the linked contig
                linked_position=contig_aln['position'],
                # First aligned base of the fragment on the linked contig
                linked_start=contig_aln['reference_start'],
                # The ID of the contig that's linked
                linked_contig=contig_aln["reference_name"],
                # The length of the linked contig
                linked_contig_len=contig_aln["reference_length"],
                # The length of the target contig
                target_contig_len=target_aln["reference_length"],
                # The ID of the read which supports this junction
                read_name=read_name,
                # Positional offset
                pos_offset=pos_offset
            )
        )

# Make a DataFrame with the results
junction_df = pd.DataFrame(junction_list)

# Write to a file
junction_df.to_csv("junctions.csv.gz", index=None)
