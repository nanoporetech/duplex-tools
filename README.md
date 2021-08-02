![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)

# Read Fillet

Read fillet is a simple utility for splitting Oxford Nanopore Sequencing reads
based on knowledge of adapter sequences. Its primary use case is for splitting
chimeric reads into their component sub-reads.


## Installation

Read Fillet can be installed from PyPI with:

    pip install read_fillet

## Usage

Read Fillet is run simply with:

    read_fillet <fastq_directory> 

To see more options run:

    read_fillet --help

For each input fastq file found in the input directory, a file with an
additional `_split` suffix will be output into the output directory, controlled
by the `--output_dir` option. If not `--output_dir` is given then the output
files are placed alongside in the input files.

The new `*_split.fastq` will contain two new reads for each read that was
split, now with suffix `_1` and `_2` Reads which were not split will also be
added to the new file, so that `*_split.fastq` can be used as a match for any
downstream analysis. For example the output may look like:

    @<read_id> <remaining_headers>
    <sequence>
    +
    <quality>
    @<read_id>_1 <remaining_headers>
    <sequence>
    +
    <quality>
    @<read_id>_2 <remaining_headers>
    <sequence>
    +
    <quality>


## Algorithmic Details

The purpose of read_fillet is to split reads likely to be concatemers of independent
reads. Such concatemers may arise from at least two mechanisms:

* chemical chimers arising from artifacts of molecule-biological steps used in the
  preparation of sequencing libraries,
* informatic chimers arising from failures of the algorithms used to process the
  primary sequencing data.

Regardless of the mechanism that has created chimeric reads, read_fillet
attempts to split such reads into their component sub-reads using knowledge of
the sequencing adapters used in library preparation. Concatemeric reads
typically contain one or more matches to a sequencing adapter internal to their
sequence; for example, commonly half-way through the read when library
read-length is from a tightly centered distribution. The following simplified
views indicate the types of error that may occur.

Here is a read, where `>>` represents the adapter sequence and `=` represent
bases from the organism of interest. The reverse-complement of the sequencing
 adapter can typically be found at the end of a read, represented as `<<` below:

```
    >>=====<<
```

A simple chimeric read containing two sub-reads can be represented in this notation as:

```
        A        B
    >>=====<<>>=====<<
```

On Oxford Nanopore sequencing platforms, and enzyme is used to unwind double-stranded
DNA and thread a single strand through a nanopore. It is possible that after the first
strand of a duplex has transited the nanopore, that the second strand is immediately
captured since it is physically close to the pore. If the capture is particularly
fast the sequencing device software does not detect a boundary between the signals
measured for the first and second strand and subsequently produced a single conjoined
read. In this case the sub-read B will be complementary to A. It is possible also that
the sub-read B is distinct from A, being derived from a second DNA duplex unrelated
to that which gave rise to read A. In this case, depending on the sequencing library,
the sequences of A and B can be unrelated.

Read-fillet looks for matches to the sequence `<<>>` and splits the original basecall
for the read into separate parts, for example the read:

Before:
```
    read_id: 0ae195a2-6993-4a0b-afa8-bb834f4739e3

            A        B
        >>=====<<>>=====<<
               ^--^
```

wil become:
```
    read_id: 0ae195a2-6993-4a0b-afa8-bb834f4739e3_1
            A
        >>=====<<

    read_id: 0ae195a2-6993-4a0b-afa8-bb834f4739e3_2
            B
        >>=====<<
```

after the splitting process. Note the `_1` and `_2` in the emitted read UUIDs.


### Assessment

An additional program `read_fillet_assess` is available to provide
an assessment of the veracity of the results provided by the main `read_fillet`
program. The assessment program uses alignment of reads to a reference sequence
to form knowledge of the true structure of reads (assuming the reference sequence
to be correct and, for example, not contain structural variants with respect to
the sequences sample).


The assessment program requires the addition dependencies: `pomoxis`,
`samtools`, and `seqkit`. These are most easily obtained using conda (or mamba
as a faster alternative):

    mamba create --name read_fillet -c bioconda seqkit samtools pomoxis python3.6
    conda activate read_fillet
    pip install read_fillet

The program examines supplementary alignments produced by the `minimap2` aligner. The
alignments are checked to determine their overlap to both the read and the reference
sequence. Reads are grouped into one of the following classes:

* **single_alignment_95%cov**: the read is associated with a single alignment spanning >95% of the reads length,
* **disjoint_with_gap**: two alignments are found for the read with an unaligned, adapter-sized gap excluded from the alignments,
* **disjoint_without_gap**: similar to the former, without an adapter-sized gap,
* **overlapping**: two alignments are found which overlap each other with respect to the reference. 
* **read_gt_2_supplementary**: the read is associated with more than two alignments.

These categories are not mutually exclusive, if a read belongs to multiple classes the priority
of labelling is from top to bottom as written above. 

To assess the veracity of the adapter-based `read_fillet` read splitting the assessment
program counts reads belonging to with- and without-gap classes as follows:

**disjoint_with_gap**

* `read_fillet` split: **true positive**
* `read_fillet` unsplit: **false negative**

**disjoint_without_gap**

* `read_fillet` split: **false positive**
* `read_fillet` unsplit: **true negative**

### Help

**Licence and Copyright**

© 2021- Oxford Nanopore Technologies Ltd.

`read_fillet` is distributed under the terms of the Mozilla Public License 2.0.

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.
