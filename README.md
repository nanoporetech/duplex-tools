![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)

# Duplex Tools

Duplex Tools contains a set of utilities for dealing with Duplex sequencing
data. Tools are provided to identify and prepare duplex pairs for basecalling
by Guppy and for recovering simplex basecalls from incorrectly concatenated
pairs.

## Installation

Duplex Tools can be installed from PyPI with:

    pip install duplex_tools

## Usage

Duplex Tools is run simply with:

    duplex_tools --help

The available sub-commands are:

* [split_on_adapter](./fillet.md) - split incorrectly concantenate duplex pairs in to their component simplex reads (formerly `read_fillet`.
* pairs_from_summary - identify candidate duplex pairs from sequencing summary output by Guppy.
* filter_pairs - filter candidate pairs using basecall-to-basecall alignment.

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
