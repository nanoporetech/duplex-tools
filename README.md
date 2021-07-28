
![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)

# Read Fillet


### Installation

Read Fillet can be installed from PyPI with:

    pip install read_fillet

### Usage

Read Fillet is run simply with:

    read_fillet <path_to_fastq_pass_directory> --pattern "*.fastq" 

To see more options run:

    read_fillet --help

For each input fastq file, a file with an additional `_split` suffix will be output into the output directory:

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
