# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
