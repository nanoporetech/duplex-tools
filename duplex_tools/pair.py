"""Pair reads.

Convenience wrapper to both form pairs (pairs_from_summary) and filter them
(filter_pairs).
"""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import pysam

import duplex_tools
from duplex_tools.filter_pairs import add_args as add_filter_args
from duplex_tools.filter_pairs import filter_candidate_pairs_by_aligning
from duplex_tools.pairs_from_summary import add_args as add_pair_args
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
                   align_threshold,
                   no_end_penalties,
                   penalty_open,
                   penalty_extend,
                   score_match,
                   score_mismatch,
                   threads,
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
    logger = duplex_tools.get_named_logger("Pair")
    find_pairs(input_bam,
               outdir=output_dir,
               max_time_between_reads=max_time_between_reads,
               max_seqlen_diff=max_seqlen_diff,
               max_abs_seqlen_diff=max_abs_seqlen_diff,
               min_qscore=min_qscore,
               )
    filter_candidate_pairs_by_aligning(f'{output_dir}/pair_ids.txt',
                                       reads_path=input_bam,
                                       bases_to_align=bases_to_align,
                                       min_length=min_length,
                                       max_length=max_length,
                                       align_threshold=align_threshold,
                                       no_end_penalties=no_end_penalties,
                                       penalty_open=penalty_open,
                                       penalty_extend=penalty_extend,
                                       score_match=score_match,
                                       score_mismatch=score_mismatch,
                                       threads=threads
                                       )

    npairs = sum(1 for _ in open(f'{output_dir}/pair_ids_filtered.txt'))
    nreads = pysam.AlignmentFile(input_bam, check_sq=False).count(
        until_eof=True)
    logger.info(f'Initial reads: {nreads}')
    logger.info(f'Created pairs: {npairs}')
    logger.info(f'Paired reads:  {2 * npairs}')
    logger.info(f'Approximate duplex rate for {input_bam}: '
                f'{2*100*npairs / nreads:.2f}%')


def argparser():
    """Create argument parser."""
    parser = ArgumentParser(
        "Filter candidate read pairs by basecall alignment.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        parents=[duplex_tools._log_level()],
        add_help=False)
    parser.add_argument(
        "bam",
        help="A bam file from dorado.")
    parser.add_argument(
        "--output_dir",
        help="The output directory", default='pairs_from_bam')

    parser = add_pair_args(parser)
    parser = add_filter_args(parser)

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
                   no_end_penalties=args.no_end_penalties,
                   align_threshold=args.align_threshold,
                   penalty_open=args.penalty_open,
                   penalty_extend=args.penalty_extend,
                   score_match=args.score_match,
                   score_mismatch=args.score_mismatch,
                   threads=args.threads,
                   )
