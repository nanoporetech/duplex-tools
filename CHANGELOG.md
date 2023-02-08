# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.3.1]
### Added
- wrapper for a prototype stereo-calling pipeline `dorado_stereo.sh` that does 2-stage stereo calling.

## [v0.3.0]
### Changed
- update fastx iteration in `split_pairs` to be compatible with pyfastx>=0.9.0.
### Fixed
- Bug where `split_pairs` would raise a StopIteration if dataset has < 5k reads. 


## [v0.2.20]
### Added
- `split_pairs`, a tool to recover non-split reads into their template/complement parts. 

## [v0.2.19]
### Fixed
- Bug where approximate start times led to incorrectly rejecting candidate duplex pairs.

## [v0.2.18]
### Fixed
- Bug where aligned BAMs would be used as inputs to filter_pairs. This meant running `filter_pairs` on a guppy directory would report an incorrect duplex rate

## [v0.2.17]
### Fixed
- Bug where an incorrect tag was used to get sequence lengths from .bam files (for dorado)
### Changed  
- ProcessPoolExecutor -> ThreadPoolExecutor in filter_pairs for faster filtering
- Reporting duplex rate in duplex_tools pair

## [v0.2.16]
### Added
- Ability to use an unmapped bam (output from dorado) as input to pairs_from_summary
- Ability to use an unmapped bam (output from dorado) as input to filter_pairs
- Convenience script (`duplex_tools pair unmapped.bam`) to call both pairs_from_summary and filter_pairs on a bam
  - usage: duplex_tools pair unmapped.bam

## [v0.2.15]
### Fixed
- Bug where template/complement pairs with small negative time between them would not be chosen (rounding error).
### Added
- Option to set --no_end_penalties in filter_pairs. This option favours partial matches and avoids unbounded negative pairing scores

## [v0.2.14]
### Fixed
- Bug where template/complement pairs with no time between them would not be chosen

## [v0.2.13]
### Added
- Option to set --threads in filter_pairs
### Updated
- Updated defaults in readme for duplex basecalling (`chunks_per_runner` 16 -> 416)

## [v0.2.12]
### Fixed
- Update defaults in pairs_from_summary to `--min_qscore 7 --max_abs_seqlen_diff 1000`

## [v0.2.11]
### Fixed
- Passed --trim_start and --trim_end from cli to main function for split_on_adapter

## [v0.2.10]
### Added
- Flags `--trim_start` and `--trim_end` for split_on_adapter (https://github.com/nanoporetech/duplex-tools/issues/13)

## [v0.2.9]
### Added
- Flag to allow splitting multiple times for reads with multiple adapters
- Moving debug output and edited reads to the main output directory

## [v0.2.8]
### Changed
- Removed explicit dependency on pathlib, which caused https://github.com/nanoporetech/duplex-tools/issues/7

## [v0.2.7]
### Fixed
- split_on_adapter also defaults to both .fastq and fastq.gz from cli
### Added
- Options --min_qscore and --max_abs_seqlen_diff to find duplex reads
- Default filtering on min_qscore (12) and maximum length difference in pairs_from_summary for duplex reads

## [v0.2.6]
### Fixed
- Surprising behaviour in split_on_adapter to only work on gzipped fastqs by default

## [v0.2.5]
### Added
- Arguments --max_length and --min_length in filter_pairs

## [v0.2.4]
### Changed
- Fixed argument bug in split_on_adapter. sample_type is now positional

## [v0.2.3]
### Changed
- Corrected documentation for fillet.md.
- Enabled multiprocessing for fastq extraction in `filter_pairs`.

## [v0.2.2]
### Changed
- Fixed number formatting in number of pairs logging.

## [v0.2.1]
### Changed
- Documentation for read splitting

## [v0.2.0]
### Changed
- Project name to duplex tools.
- Create single entry point program with subcommands.
### Added
- Duplex read pairing and filtering programs.

## [v0.1.3]
### Added
- Integration test to run the main read_fillet entry point.
### Fixed
- Bug that caused compressed outputs to be missing the .gz extension.

## [v0.1.2]
### Fixed
- Bug that caused setting trim=0 to return an empty sequence.

## [v0.1.1]
### Changed
- Updated README.

## [v0.1.0]
### Added
- First version.
