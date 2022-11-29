"""Pair reads.

Convenience wrapper to both form pairs (pairs_from_summary) and filter them
(filter_pairs).
"""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import duplex_tools
from duplex_tools.filter_pairs import filter_candidate_pairs_by_aligning
from duplex_tools.pairs_from_summary import find_pairs


def pair_and_align(input_bam,
                   max_time_between_reads,
                   max_seqlen_diff,
                   max_abs_seqlen_diff,
                   min_qscore,
                   bases_to_align,
                   min_length,
                   max_length,
                   output_dir,
                   **kwargs):
    """Pair and align reads from an unmapped bam.

    :param input_bam: The input bam file (containing unmapped reads)
    :param output_dir: The output directory (the pair_ids_filtered.txt is here)
    :param max_time_between_reads: see pairs_from_summary
    :param max_seqlen_diff: see pairs_from_summary
    :param max_abs_seqlen_diff: see pairs_from_summary
    :param min_qscore: see pairs_from_summary
    :param bases_to_align: see filter_pairs
    :param min_length: see filter_pairs
    :param max_length: see filter_pairs
    """
    reads_directory = str(Path(input_bam).parent)
    find_pairs(input_bam,
               outdir=output_dir,
               max_time_between_reads=max_time_between_reads,
               max_seqlen_diff=max_seqlen_diff,
               max_abs_seqlen_diff=max_abs_seqlen_diff,
               min_qscore=min_qscore,
               )
    filter_candidate_pairs_by_aligning(f'{output_dir}/pair_ids.txt',
                                       reads_directory,
                                       bases_to_align=bases_to_align,
                                       min_length=min_length,
                                       max_length=max_length,
                                       )


def argparser():
    """Create argument parser."""
    parser = ArgumentParser(
        "Filter candidate read pairs by basecall alignment.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        parents=[duplex_tools._log_level()], add_help=False)
    parser.add_argument(
        "bam",
        help="A bam file from dorado.")
    parser.add_argument(
        "--output_dir",
        help="The output directory", default='pairs_from_bam')
    parser.add_argument(
        "--max_time_between_reads", type=int, default=10000,
        help=(
            "Maximum time (seconds) between reads for them to be "
            "deemed a pair."))
    parser.add_argument(
        "--max_seqlen_diff", type=float, default=0.4,
        help=(
            "Maximum ratio (a - b) / a, where a and b are the "
            "sequence lengths of a putative pair."))
    parser.add_argument(
        "--max_abs_seqlen_diff", type=int, default=50000,
        help=(
            "Maximum sequence length difference between template and "
            "complement"))
    parser.add_argument(
        "--min_qscore", type=float, default=0,
        help=(
            "The minimum simplex qscore required from both template and "
            "complement"))
    parser.add_argument(
        "--bases_to_align", default=250, type=int,
        help="Number of bases from each read to attempt alignment.")
    parser.add_argument(
        "--min_length", default=1, type=int,
        help="Minimum length of template and complement.")
    parser.add_argument(
        "--max_length", default=float('inf'), type=float,
        help="Maximum length of template and complement.")
    return parser


def main(args):
    """Entry point."""
    pair_and_align(input_bam=args.bam,
                   output_dir=args.output_dir,
                   max_time_between_reads=args.max_time_between_reads,
                   max_seqlen_diff=args.max_seqlen_diff,
                   max_abs_seqlen_diff=args.max_abs_seqlen_diff,
                   min_qscore=args.min_qscore,
                   bases_to_align=args.bases_to_align,
                   min_length=args.min_length,
                   max_length=args.max_length,
                   )
