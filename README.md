
![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)

# Read Fillet


### Installation

Read Fillet can be installed from PyPI with:

    pip install read_fillet

### Usage

Read Fillet is run simply with:

    read_fillet <fastq_directory> 

To see more options run:

    read_fillet --help

For each input fastq file found in the input directory, a file with an additional `_split` suffix will be output into the output directory, controlled by the `--output_dir` option. If not `--output_dir` is given then the output files are placed alongside in the input files.

The new `*_split.fastq` will contain two new reads for each read that was split, now with suffix `_1` and `_2`
Reads which were not split will also be added to the new file, so that `*_split.fastq` can be used as a match for any downstream analysis. For example the output may look like:

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

### Extra Details


The purpose of this script is to split reads that are likely to be informatic
concatamers.  These reads typically have an imperfect match to two neighbouring
adapter sequences found next to each other.

A simplified view is the following:

Here is a read, where ">>" represents the adapter sequence and "=" represent
bases from the organism of interest.  The reverse adapter can typically be
found in the end of a read, represented as "<<" below

A typical read, with adapters in the ends, and bases of interest in the middle

```
    >>=====<<
```

Here is a chimeric read that has not been split in the software.

```
        A        B
    >>=====<<>>=====<<
```

B can either be complementary to A, in which case they form a follow-on pair,
which formed dsDNA prior to sequencing.  B can also be genomically distant from
A, in which case it is may also desirable to split this read to avoid chimaeras

This module will look for matches to "<<>>" like below, and split these reads
out into separate parts.

Before:
```
    read_id: 0ae195a2-6993-4a0b-afa8-bb834f4739e3

            A        B
        >>=====<<>>=====<<
               ^--^
```

After:
```
    read_id: 0ae195a2-6993-4a0b-afa8-bb834f4739e3_1
            A
        >>=====<<

    read_id: 0ae195a2-6993-4a0b-afa8-bb834f4739e3_2
            B
        >>=====<<
```

**Assessment**

An additional assessment program `read_fillet_assess` is available. This requires the addition dependencies: `pomoxis`, `samtools`, and `seqkit`. These are most easily obtained using conda (or mamba as a faster alternative):

    mamba create --name read_fillet -c bioconda seqkit samtools pomoxis python3.6
    conda activate read_fillet
    pip install read_fillet


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
