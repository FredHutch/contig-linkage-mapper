# Contig Linkage Mapper

#### Quantify the rate of alternate contig linkages from whole-genome shotgun data by alignment

## Background

One of the challenges of analyzing _de novo_ assemblies of complex genomes is that mobile or repeated
genetic elements will often break assembly algorithms. In the context of metagenomic analysis, this
phenomenon is often observed for transposons, bacteriophages, and other elements which may be present
between different genomic flanking regions when comparing different microbial strains or isolates.
From the perspective of information content, mobile and repeated genetic elements will always present
a challenge for genomic and metagenomic analysis in any situation where the size of the genetic element
is longer than the length of the sequencing reads which can be generated.

While long-read sequencing is becoming increasingly accurate and affordable, short-read sequencing
platforms are extremely robust and mature. With long reads (e.g. >1kb), mobile elements can be linked
to flanking regions by identifying reads which span the junction with the adjacent genomic context.
With short reads, such linkage information is generated when the two ends of a sequenced molecule
are aligned on either side of the junction.

## Approach

The use-case which is targeted by this utility is when the user wishes to identify the amount of
evidence supporting the linkage of a single contig of interest to different adjacent genomic elements
across different biological specimens. The input data expected by this utility consists of:

1. Paired-end short-read whole-genome shotgun sequencing (WGS) datasets in FASTQ format
2. The DNA sequence of a single genetic element of interest in FASTA format
3. The DNA sequences of _de novo_ assembled contigs from the same specimens in FASTA format

Using this set of input data, this utility will:

a. Align WGS data (1) against both the target contig (2) and the larger assemblies (3)
b. Identify reads for which only one half of the pair is aligned to the contig of interest
c. Summarize the placement of the reads paired to those identified in (b) on contigs from the larger assembly (3)

### Inferring WGS insert size

The breakpoint of an assembly can have unresolved spacing of adjacent sequences. However, paired-end
sequencing data is often generated with a constrained distribution of DNA molecule lengths. That
distribution can be inferred from read alignment datasets based on the relative position of those
reads for which both ends align in the expected orientation to the same contig. By applying that
"insert size" data to the contig junction mapping, we can infer the general offset of a contig
junction which is supported by the available sequencing data.

### Summarizing junctions

Each pair of reads which provides evidence linking the target contig to another is considered
to be a 'junction', or a fragment of DNA which spans a contig-contig junction. While there
are many different arrangements and orientations which may be observed, we will attempt to
simplify the reporting of these data by using the concepts of 'internal' and 'end' junctions.

In this analysis, an 'end' junction is one which is estimated to span either the 5' or 3'
terminus of the contig. To determine the distance threshold which a read alignment must fall
within in order to be classified as an 'end' alignment, we will use the distribution of insert
sizes calculated from the paired-end alignments to the target contig (described above).
The default behavior for the analysis will be mark a read alignment as 'end' if it is no more
than `N` bases from the end of the contig which it is oriented towards (3' for forward alignments
and 5' for reverse alignments), where `N` is the mean insert size + `sigma` standard deviations.
The parameter `sigma` is set by default at 1, but may be modified by the user.

The labels of 'internal' and 'end' will be applied to both the target contig and the linked contig.

Junctions will be summarized from the perspective of the read which aligns to the target contig.
Reads which align in the forward orientation are pointing towards the 3' end of the contig, with
reverse strand alignment pointing towards the 5' end of the contig. If the paired read on the
linked contig is aligned to the opposite strand, then the contig orientation will be marked as
'fwd'. If the pair of reads are aligned on the same strand of the two contigs, the contig
orientation will be marked as 'rev'.

Finally, the `offset` of the junction will indicate the position of the linked contig relative
to the target contig, applying the median insert size to the paired read positions. For 'fwd'
orientations, the offset will apply to the 5' end of the paired contig, while 'rev' orientations
will use the 3' end of the contig. The motivation for using this type of coordinate system is
that any linked contig which is expected to directly abut the target contig (based on the position
of the aligned reads) will have an `offset` value which is either very close to the length of the
target contig (for 'fwd' orientations) or the negative value of the length of the linked contig
(for 'rev' orientations).

For a visual demonstration of the contig offset concept, refer to the diagrams included in
`docs/positional.offset.pttx`.

## Output

The output generated by this utility will be summarized in terms of:

- The direction of the linkage, either from the 'fwd' or 'rev' strand of the target contig (2)
- The position of the junction in the target contig, 'internal' or 'end'
- The 'linked' contig from the larger assembly (3)
- The position of the junction in the linked contig, 'internal' or 'end'
- The orientation of linked contig relative to the target contig ('fwd' or 'rev')
- The length of the linked contig
- The linked contig positional offset
- The name of the specimen from which reads were detected supporting this junction
- The number of reads from that dataset which support this junction

NOTE: Reads which are consistent with a circular structure for the contig of
interest will be included in the output above, using the ID `self` for the linked contig.
The total size of the circular moledule can be inferred by applying the contig junction
offset to the contig length.

Other outputs which will be generated to help put the data in context are:

- The total number of reads from each specimen
- The number of reads from each specimen which aligned to the contig of interest

## Analysis Details

The workflow will take the approach of:

1. Align input reads to the contig of interest
2. Segment the aligned reads into those which:
- (2a) Align both ends in proper orientation
- (2b) Align both ends in outward-facing orientation (i.e., circular orientation)
- (2c) Align only one end to the target contig
3. Infer the per-specimen DNA insert size from (2a)
4. Align the other pair of the reads from (2c) to the larger set of _de novo_ assembled contigs (input 3, above)
5. Aggregate all contig linkage information and generate outputs
