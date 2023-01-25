![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)

# Duplex Tools

Duplex Tools contains a set of utilities for dealing with Duplex sequencing
data. Tools are provided to identify and prepare duplex pairs for basecalling
by Dorado (recommended) and Guppy, and for recovering simplex basecalls from incorrectly concatenated
pairs.

## Installation

Duplex Tools is written in Python and can be installed directly from PyPI.
We recommend installing Duplex Tools into an isolated virtual environment
by following:

    python -m venv venv --prompt duplex
    . venv/bin/activate
    pip install duplex_tools

after which the code tools will be available using the `duplex_tools` command.

## General Usage

Duplex Tools is run simply with:

    duplex_tools --help

The available sub-commands are:

### Duplex pairing


#### Compatible with Dorado
* `pair` - a wrapper to pair duplex reads, using `pairs_from_summary` and then `filter_pairs`.
* `split_pairs` - a utility for recovering and pairing duplex reads (for cases where template/complement are contained within a single minknow read).

#### Compatible with Guppy+Dorado
* `pairs_from_summary` - identify candidate duplex pairs from sequencing summary output by Guppy or unmapped SAM/BAM by dorado.
* `filter_pairs` - filter candidate pairs using basecall-to-basecall alignment.

### Additional tools
* [split_on_adapter](./fillet.md) - split the non-split duplex pairs in to their component simplex reads (formerly `read_fillet`). 
  * This tool splits basecalled sequences into new sequences. For this reason, it's possible to perform _basespace_ duplex calling after using this method, but not regular stereo calling


## Usage with Dorado (recommended)

Currently, pairing and calling are separate processes to allow for workflow flexibility.

For greatest duplex recovery, follow these steps:

1) Simplex basecall with dorado (with `--emit-moves`)
2) Pair reads
3) Duplex-basecall reads


### 1a) Simplex basecall with dorado
This will create an (unmapped) .sam file which has a mapping between the signal and bases.
`--emit-moves` allows for additional pairs to be found in step 2b.

    $ dorado basecaller dna_r10.4.1_e8.2_400bps_fast@v4.0.0 pod5s/ --emit-moves > unmapped_reads_with_moves.sam

### 2a) Find duplex pairs for Dorado stereo/basespace basecalling
This will detect the majority of pairs and put them in the `pairs_from_bam` directory.

    duplex_tools pair --output_dir pairs_from_bam unmapped_reads_with_moves.bam


### 2b) Find additional duplex pairs in non-split reads (optional)

The steps below can recover non-split pairs and allows duplex-calling of them.

**Use the sam and a pod5 directory to create additional pairs**

    $ duplex_tools split_pairs unmapped_reads_with_moves.sam pod5s/ pod5s_splitduplex/
    $ cat pod5s_splitduplex/*_pair_ids.txt > split_duplex_pair_ids.txt

### 3) Stereo basecall all the reads

From the main pairing:

    $ dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.0.0 pod5s/ --pairs pairs_from_bam/pair_ids_filtered.txt > duplex_orig.sam

From the additional pairing (optional):

    $ dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.0.0 pod5s_splitduplex/ --pairs split_duplex_pair_ids.txt > duplex_splitduplex.sam


## Usage with Guppy

**Preparing duplex reads for Guppy duplex basecalling**

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
        -x 'cuda:0' -c dna_r10.4.1_e8.2_400bps_sup.cfg \
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
