![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)

# Duplex Tools

Duplex Tools contains a set of utilities for dealing with Duplex sequencing
data. Tools are provided to identify and prepare duplex pairs for basecalling
by Guppy and for recovering simplex basecalls from incorrectly concatenated
pairs.

## Installation

Duplex Tools is written in Python and can be installed directly from PyPI.
We recommend installing Duplex Tools into an isolated virtual environment
by following:

    python -m venv venv --prompt duplex
    . venv/bin/activate
    pip install duplex_tools

after which the code tools will be available using the `duplex_tools` command.

## Usage

Duplex Tools is run simply with:

    duplex_tools --help

The available sub-commands are:

* [split_on_adapter](./fillet.md) - split incorrectly concantenate duplex pairs in to their component simplex reads (formerly `read_fillet`).
* pairs_from_summary - identify candidate duplex pairs from sequencing summary output by Guppy.
* filter_pairs - filter candidate pairs using basecall-to-basecall alignment.
* pair - a wrapper to `pairs_from_summary` and then `filter_pairs`. Currently only compatible with dorado 


**Preparing duplex reads for Guppy basecalling**

To prepare reads for duplex calling Duplex Tools provides two programs. The
first parses the sequencing summary output by the Guppy basecaller (or the metadata in a .bam or .sam from dorado) in order
to generate candidate pairs from examining simple read metrics. The second
program analyses the basecalls of candidate reads, checking for similarity.

To run the basic sequencing summary(/bam metadata) based pairing run the following:

    duplex_tools pairs_from_summary <sequencing_summary.txt/dorado.bam> <output directory>

The primary output of the above will be a text file named `pair_ids.txt` in the
user specified output directory. Although this file can be given to Guppy to perform
duplex calling we recommend running the second basecall-to-basecall alignment
filtering provided by the `filter_pairs` command:

    duplex_tools filter_pairs <pair_ids.txt> <fastq directory/dorado.bam>

The first option here is the file described above and output by `pairs_from_summary`.
The second option should be specified as the Guppy (or MinKNOW), or dorado output directory
containing `fastq` or `bam` data --- the directory will be search recursively for all `.fastq.gz`, `.fastq`, and `.sam/.bam` files. 

The output of this second command will be a file named
`pair_ids_filtered.txt` placed alongside the `pair_ids.txt` file.

**Duplex basecalling with Guppy**

The file `pair_ids_filtered.txt` as prepared above can be used with the
original `.fast5`/`.pod5` files produced during a sequencing run in order to calculate
high quality duplex basecalls.

For example,

    guppy_basecaller_duplex \
        -i <MinKNOW directory> \
        -r -s duplex_calls \
        -x 'cuda:0' -c q20-fixed-2.0-ft-10M.cfg \
        --chunks_per_runner 416 \
        --duplex_pairing_mode from_pair_list \
        --duplex_pairing_file pair_ids_filtered.txt

will produce duplex basecalls using the read pairs stored in the
`pair_ids_filtered.txt` file using `.fast5`/`.pod5` files found in the user
provided MinKNOW output directory.

**Duplex basecalling with Dorado**

Please use `duplex_tools pair unmapped_dorado.bam`. 
This will run both the pairing and pairwise alignment-based filtering to get a pair_ids_filtered.txt that can be passed to dorado. 


For more details, see https://github.com/nanoporetech/dorado. 


### Help

**Licence and Copyright**

© 2021- Oxford Nanopore Technologies Ltd.

`Duplex Tools` is distributed under the terms of the Mozilla Public License 2.0.

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
