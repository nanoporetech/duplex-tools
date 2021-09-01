"""Filter candidate pairs by inspecting mutual basecall alignment."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from collections import defaultdict
import glob
from pathlib import Path
import pickle

import pandas as pd
import parasail
import pysam

import duplex_tools


comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    # '-': '-'
}
comp_trans = str.maketrans(''.join(comp.keys()), ''.join(comp.values()))


def reverse_complement(seq):
    """Reverse complement sequence.

    :param: input sequence string.

    :returns: reverse-complemented string.
    """
    return seq.translate(comp_trans)[::-1]


def filter_candidate_pairs_by_aligning(
        read_pairs: str,
        fastq: str,
        bases_to_align: int = 250,
        align_threshold: float = 0.6,
        penalty_open: int = 4,
        penalty_extend: int = 1,
        score_match: int = 2,
        score_mismatch: int = -1,
        loglevel: str = "INFO") -> None:
    """Filter candidate read pairs by quality of alignment.

    :param read_pairs: Path to file with two space-separated read-ids per row,
        the leftmost coming first in time.
    :param fastq: The path to a fastq file with _all_ reads, both passing and
        failing.
    :param bases_to_align: Number of bases to use from end of the first read,
        and from the start of second.
    :param align_threshold: Which alignment threshold to use for passing
        pairs.
    :param penalty_open: Open penalty passed to parasail.
    :param penalty_extend: Extend penalty passed to parasail.
    :param score_match: Match score passed to parasail.
    :param score_mismatch: Mismatch score passed to parasail.
    :param loglevel: Level to log at, for example 'INFO' or 'DEBUG'.

    This function takes a path to a file with pairs of candidate followon
    read-ids and require a fastq that contains the same read-ids.

    The end of the first read is aligned to the reverse complement of the
    beginning of the second read, as shown below:
    """
    logger = duplex_tools.get_named_logger("FiltPairs")
    logger.info(
        f"Filtering Parameters: "
        f"\n\tnbases_per_strand_to_align={bases_to_align}"
        f"\n\talign_threshold={align_threshold}")
    logger.info(
        f"Alignment Parameters: "
        f"\n\tscore_match={score_match}"
        f"\n\tscore_mismatch={score_mismatch}"
        f"\n\tpenalty_open={penalty_open}"
        f"\n\tpenalty_extend:{penalty_extend}")
    read_pairs = Path(read_pairs)
    # Index and read pairs
    pairs = pd.read_csv(read_pairs, sep=" ", names=["first", "second"])
    if fastq.endswith(".pkl"):  # for testing
        logger.info("Extracting read end data from pickle file.")
        with open(fastq, 'rb') as fh:
            fastq_index = pickle.load(fh)
    else:
        fastq_index = read_all_fastq(
            fastq, pairs, bases_to_align)
        # dump to pickle
        pkl = Path(read_pairs.parent, "read_segments.pkl")
        with open(pkl, "wb") as fh:
            pickle.dump(fastq_index, fh)
    logger.info("Starting alignments.")

    # Align all of them
    alignment_scores_df = align_all_pairs(
        align_threshold, fastq_index, bases_to_align, pairs,
        penalty_extend, penalty_open, score_match, score_mismatch)

    # Finally, write full summary and filtered pairs
    alignment_scores_df.to_csv(
        Path(read_pairs.parent, f"{read_pairs.stem}_scored.csv"),
        index=False)
    alignment_pairs_filtered = alignment_scores_df.query(
        f"score > {align_threshold}")
    alignment_pairs_filtered[["read_id", "read_id_next"]].to_csv(
        Path(read_pairs.parent, f"{read_pairs.stem}_filtered.txt"),
        index=False, header=False, sep=" ")


def read_all_fastq(fastq, pairs, n_bases):
    """Find an read all necessary data from fastq files."""
    logger = duplex_tools.get_named_logger("ReadFastq")
    first = set(pairs["first"])
    second = set(pairs["second"])

    def _get_files():
        for ext in ("fastq", "fastq.gz", "fq", "fq.gz"):
            files = glob.iglob(
                "{}/**/*.{}".format(fastq, ext), recursive=True)
            yield from files

    results = dict()
    files = list(_get_files())
    for i, file in enumerate(files):
        if i % 50 == 0:
            logger.info(
                "Processed {}/{} input fastq files.".format(i, len(files)))
        logger.debug("Extracting read ends from: {}".format(file))
        for read in pysam.FastxFile(file, persist=False):
            if read.name in first:
                results[(read.name, 0)] = str(read.sequence[-n_bases:])
            if read.name in second:  # a read can be in both
                results[(read.name, 1)] = reverse_complement(
                    str(read.sequence[:n_bases]))
    complete = 100 * len(results) / (2 * len(pairs))
    logger.info(
        "Found {:.1f}% of required reads.".format(complete))
    return results


def align_all_pairs(
        align_threshold,
        fastq_index,
        bases_to_align,
        pairs,
        penalty_extend,
        penalty_open,
        score_match,
        score_mismatch) -> pd.DataFrame:
    """Align read pairs to each other using parasail."""
    counter = defaultdict(int)
    alignment_scores = list()
    npairs = len(pairs)
    logger = duplex_tools.get_named_logger("AlignPairs")
    logger.info(f"Aligning {npairs} pairs")
    score_matrix = parasail.matrix_create("ACGT", score_match, score_mismatch)
    for idx, read_pair in enumerate(pairs.itertuples(), 1):
        if idx % 50000 == 0:
            done = 100 * idx / npairs
            good = 100 * counter["good"] / idx
            skip = 100 * counter["skipped"] / idx
            logger.info(
                f"Processed/Skip/Good: {done:.0f}%/{skip:.0f}%/{good:.0f}%")

        try:
            seq1 = fastq_index[(read_pair.first, 0)]
        except KeyError:
            logger.debug(
                f"Skipped {read_pair.first}: sequence missing.")
            counter["skipped"] += 1
            counter["read0 missing"] += 1
            continue
        try:
            seq2 = fastq_index[(read_pair.second, 1)]
        except KeyError:
            logger.debug(
                f"Skipped {read_pair.second}: sequence missing.")
            counter["skipped"] += 1
            counter["read1 missing"] += 1
            continue

        # TODO: why is this necessary, move to read_all_fastq
        if len(seq1) == 0 or len(seq2) == 0:
            logger.debug(f"Skipped {read_pair}, reads too short.")
            counter["skipped"] += 1
            continue

        # Run a semi-global alignment (sg) with zero end-penalty for seq2 (dx)
        result = parasail.sg_dx_trace_scan_16(
            seq1, seq2, penalty_open, penalty_extend, score_matrix)

        # scale score to read length
        score_followon = result.score / result.len_ref
        if score_followon > align_threshold:
            logger.debug(
                f"\n{read_pair.first} {read_pair.second}: {score_followon}")
            logger.debug(f"\n{seq1}\n{seq2}")

        alignment_scores.append(
            (read_pair.first, read_pair.second, score_followon))
        align_quality = "good" if score_followon > align_threshold else "bad"
        counter[align_quality] += 1

    alignment_scores_df = pd.DataFrame(
        alignment_scores, columns=["read_id", "read_id_next", "score"])
    logger.info("Good pairs: {}".format(counter["good"]))
    logger.debug(counter)
    return alignment_scores_df


def argparser():
    """Create argument parser."""
    parser = ArgumentParser(
        "Filter candidate read pairs by basecall alignment.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        parents=[duplex_tools._log_level()], add_help=False)
    parser.add_argument(
        "read_pairs",
        help="Candidate space-separated read ID pairs, time ordered.")
    parser.add_argument(
        "fastq",
        help="Directory to search of fastq(.gz) files.")
    parser.add_argument(
        "--bases_to_align", default=250, type=int,
        help="Number of bases from each read to attempt alignment.")
    parser.add_argument(
        "--loglevel",
        help="Level to log at, for example 'INFO' or 'DEBUG'.")
    grp = parser.add_argument_group("score options")
    grp.add_argument(
        "--align_threshold", default=0.6, type=float,
        help="Alignment score threshold (per-base) for pairing decision.")
    grp.add_argument(
        "--penalty_open", default=4, type=int,
        help="Open penalty passed to parasail.")
    grp.add_argument(
        "--penalty_extend", default=1, type=int,
        help="Extend penalty passed to parasail.")
    grp.add_argument(
        "--score_match", default=2, type=int,
        help="Match score passed to parasail.")
    grp.add_argument(
        "--score_mismatch", default=-1, type=int,
        help="Mismatch score passed to parasail.")
    return parser


def main(args):
    """Entry point."""
    filter_candidate_pairs_by_aligning(
        args.read_pairs, args.fastq,
        args.bases_to_align, args.align_threshold,
        args.penalty_open, args.penalty_extend,
        args.score_match, args.score_mismatch)
