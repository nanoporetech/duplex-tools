"""Split template/complement pairs which have come out as a single read."""
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import os

import duplex_tools
from duplex_tools.split_pairs_steps import get_split_points, split_pod5


def split_pairs(
    alignment_file: str,
    pod5_input_dir: str,
    pod5_output_dir: str,
    num_threads: int,
    max_reads: int = None,
    force_overwrite=False,
    debug_dir=None,
    match_threshold=0.8,
    left_midpoint_threshold=0.45,
    right_midpoint_threshold=0.55
):
    """Detect template/complement pairs which have come out as single reads.

    :param left_midpoint_threshold: Require the midpoint to be -> of here.
    :param right_midpoint_threshold: Require the midpoint to be <- of here.
    :param alignment_file: The uBAM/uSAM to use for splitting
                           (needs moves in mv tag)
    :param pod5_input_dir: The input directory
    :param pod5_output_dir: The output directory
    :param num_threads: The number of threads to use
    :param max_reads: Maximum number of reads to write out (for debug)
    :param force_overwrite: Whether to force overwrite of new pod5 files
    :return:
    """
    # First detect the locations where the reads should be split
    # The "split_locations" will be a dictionary like:
    # [{'read_id': {'left':(1,10), 'right': (20,30)}}]
    logger = duplex_tools.get_named_logger("SplitPairs")
    logger.info('Find split locations.')
    logger.info(f'Self-align reads from: {alignment_file}')
    split_locations = get_split_points(
        alignment_file,
        num_threads,
        max_reads=max_reads,
        match_threshold=match_threshold,
        left_midpoint_threshold=left_midpoint_threshold,
        right_midpoint_threshold=right_midpoint_threshold
    )

    # Secondly, use the split locations to extract the reads into new pod5
    # files
    logger.info(f'Splitting {len(split_locations)} reads from: '
                f'{pod5_input_dir} into {pod5_output_dir}')
    _ = split_pod5(
        pod5_input_dir,
        split_locations,
        pod5_output_dir,
        force_overwrite=force_overwrite,
        debug_dir=debug_dir
    )


def argparser():
    """Create argument parser."""
    parser = ArgumentParser(
        "Filter candidate read pairs by basecall alignment.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        parents=[duplex_tools._log_level()],
        add_help=False,
    )
    parser.add_argument(
        "alignment_file",
        help="A uBAM or uSAM file from dorado containing moves ",
    )
    parser.add_argument(
        "pod5_input_dir", help="A directory containing at least one pod5 file."
    )
    parser.add_argument(
        "pod5_output_dir",
        help="A directory where the new (split) reads will be placed",
    )
    parser.add_argument(
        "--max_reads",
        type=int,
        default=None,
        help="The maximum number of reads to process",
    )
    parser.add_argument(
        "--force_overwrite",
        action="store_true",
        help="Whether to force overwriting existing pod5s",
    )
    parser.add_argument(
        "--debug_dir",
        type=str,
        default=None,
        help="Which directory to write plots of raw signal to",
    )
    parser.add_argument(
        "--match_threshold",
        type=float,
        default=0.6,
        help="Require this fraction of the template to align to the "
             "complement",
    )
    parser.add_argument(
        "--left_midpoint_threshold",
        type=float,
        default=0.45,
        help="Require this fraction of the template to align to the "
             "complement",
    )
    parser.add_argument(
        "--right_midpoint_threshold",
        type=float,
        default=0.55,
        help="Require this fraction of the template to align to the "
             "complement",
    )
    ncpu = os.cpu_count()
    if ncpu > 2:
        nthreads_default = ncpu-1
    else:
        nthreads_default = ncpu
    parser.add_argument("--threads",
                        type=int,
                        help="Number of threads to use",
                        default=nthreads_default)

    return parser


def main(args):
    """Entry point."""
    split_pairs(
        alignment_file=args.alignment_file,
        pod5_input_dir=args.pod5_input_dir,
        pod5_output_dir=args.pod5_output_dir,
        num_threads=args.threads,
        max_reads=args.max_reads,
        force_overwrite=args.force_overwrite,
        debug_dir=args.debug_dir,
        match_threshold=args.match_threshold,
        left_midpoint_threshold=args.left_midpoint_threshold,
        right_midpoint_threshold=args.right_midpoint_threshold
    )
