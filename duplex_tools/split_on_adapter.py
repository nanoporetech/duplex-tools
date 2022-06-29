"""Split reads containing internal adapter sequences."""
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from concurrent.futures import ProcessPoolExecutor
import functools
import gzip
from pathlib import Path
import pickle
import sys

import edlib
from more_itertools import pairwise
from natsort import natsorted
import numpy as np
from pyfastx import Fastx
from tqdm import tqdm

import duplex_tools

EDIT_THRESHOLDS = {'PCR': 45, 'Native': 9}
mask_size_default_head = 5
mask_size_default_tail = 14
mask_size_default_N = 11
HEAD_ADAPTER = 'AATGTACTTCGTTCAGTTACGTATTGCT'
TAIL_ADAPTER = 'GCAATACGTAACTGAACGAAGT'
rctrans = str.maketrans('ACGT', 'TGCA')


def rev_comp(seq):
    """Reverse complement a DNA sequence."""
    return str.translate(seq, rctrans)[::-1]


def build_targets(
        n_bases_to_mask_head, n_bases_to_mask_tail, degenerate_bases=None,
        pcr_primers=(
            'ACTTGCCTGTCGCTCTATCTTCGGCGTCTGCTTGGGTGTTTAACC',
            'TTTCTGTTGGTGCTGATATTGCGGCGTCTGCTTGGGTGTTTAACCT'),
        head_adapter=HEAD_ADAPTER, tail_adapter=TAIL_ADAPTER,
        n_replacement=None):
    """Build some targets."""
    if degenerate_bases is None:
        degenerate_bases = n_bases_to_mask_head + n_bases_to_mask_tail

    if n_replacement is None:
        middle_sequence = 'N' * degenerate_bases
    else:
        middle_sequence = n_replacement

    targets = {
        'Native': [
            (tail_adapter[:len(tail_adapter) - n_bases_to_mask_tail]
                + middle_sequence
                + head_adapter[n_bases_to_mask_head:])],
        'PCR': [
            (rev_comp(x)
                + tail_adapter[:len(tail_adapter) - n_bases_to_mask_tail]
                + middle_sequence + head_adapter[n_bases_to_mask_head:] + y)
            for x in pcr_primers for y in pcr_primers]
    }
    return targets


def find_mid_adaptor(
        seq, targets, print_alignment=False, print_threshold=10,
        print_id=None, trim_start=200, trim_end=200):
    """Find adapters in middle of reads."""
    seq = seq[trim_start:-trim_end or None]  # remove start and end adaptor
    results = [
        edlib.align(
            target, seq, mode="HW", task="path",
            additionalEqualities=(
                ('N', 'A'), ('N', 'C'), ('N', 'G'), ('N', 'T')))
        for target in targets]
    if print_alignment:
        alignments = [
            edlib.getNiceAlignment(result, target, seq)
            for result, target in zip(results, targets)
            if result['cigar'] is not None]
        for alignment, result in zip(alignments, results):
            if result['editDistance'] < print_threshold:
                print(f"{print_id} editdistance-{result['editDistance']}")
                print("\n".join(alignment.values()))

    i = np.argmin([x['editDistance'] for x in results])
    res = results[i]
    if res['cigar'] is not None:
        res['locations'] = [
            (x + trim_start, y + trim_start)
            for (x, y) in res['locations']]
    return res


def write_match_to_fasta(file, seq, start, end, read_id):
    """Write matches to fasta."""
    matchlen = end - start
    matchlen_half = matchlen // 2
    seq_left = seq[(start + matchlen_half - 100):start]
    seq_mid = seq[start:end]
    seq_right = seq[end:(end - matchlen_half + 100)]
    fullseq = ''.join([seq_left, seq_mid, seq_right])
    if len(fullseq) > 0:
        file.write(f'>{read_id}\n{fullseq}\n')


def deduplicate_locations_first_key(result):
    """Deduplicate locations."""
    result['locations'] = sorted(
        list({x: y for x, y in reversed(result['locations'])}.items()))
    return result


def process_file(
        fastx, targets, output_dir=None,
        debug_output=False,
        edit_threshold=None,
        print_alignment=False,
        print_threshold_delta=0,
        allow_multiple_splits=False,
        trim_start=200,
        trim_end=200
        ):
    """Run the workflow on a single file."""
    newfastx = fastx.with_name(
        fastx.name.replace('.fastq',
                           '').replace('.gz',
                                       '') + '_split').with_suffix('.fastq.gz')
    if output_dir is not None:
        newfastx = Path(output_dir) / newfastx.name
    if debug_output:
        newfasta = Path(output_dir) / Path(
            fastx.stem.split('.')[0] + '_middle').with_suffix('.fasta')
        fasta = open(newfasta, 'w')
    edited_reads = set()
    unedited_reads = set()
    split_multiple_times = set()

    nwritten = 0
    with gzip.open(newfastx, mode='wt', compresslevel=1) as outfh:

        for read_id, seq, qual, comments in \
                tqdm(Fastx(str(fastx)), leave=False):
            result = find_mid_adaptor(
                seq, targets,
                print_alignment=print_alignment,
                print_threshold=edit_threshold + print_threshold_delta,
                print_id=read_id,
                trim_start=trim_start,
                trim_end=trim_end)

            if result['editDistance'] < edit_threshold:
                result = deduplicate_locations_first_key(result)
                if not allow_multiple_splits and len(result['locations']) > 1:
                    outfh.write(f'@{read_id} {comments}\n{seq}\n+\n{qual}\n')
                    split_multiple_times.add(read_id)
                    nwritten += 1
                    continue
                else:
                    hits = []
                    edited_reads.add(read_id)
                    for left_hit, right_hit in pairwise(
                            [(0, 0), *result['locations'], (len(seq),
                                                            len(seq))]):
                        hits.append([left_hit[1], right_hit[0]])
                    for idx, (start, end) in enumerate(hits, start=1):
                        if debug_output:
                            write_match_to_fasta(fasta,
                                                 seq,
                                                 start,
                                                 end,
                                                 read_id)
                        # This edge case can happen and results in empty
                        # sequence
                        if end <= start:
                            continue
                        subseq = seq[start:end]
                        subqual = qual[start:end]
                        h = (f'@{read_id}_{idx} {comments} {start}->{end}\n'
                             f'{subseq}\n'
                             f'+\n'
                             f'{subqual}\n')
                        outfh.write(h)
                        nwritten += 1
            else:
                outfh.write(f'@{read_id} {comments}\n{seq}\n+\n{qual}\n')
                unedited_reads.add(read_id)
    if debug_output:
        fasta.close()
    return edited_reads, unedited_reads, split_multiple_times, nwritten


def split(
        fastq_dir,
        output_dir=None,
        type="Native",
        n_bases_to_mask_tail=mask_size_default_tail,
        n_bases_to_mask_head=mask_size_default_head,
        degenerate_bases=mask_size_default_N,
        debug_output=False,
        edit_threshold=None,
        n_replacement=None,
        pattern='*.fastq*', print_alignment=False,
        print_threshold_delta=0,
        threads=None,
        allow_multiple_splits=False,
        trim_start=200,
        trim_end=200,
        ):
    """Split reads.

    :param fastq_dir: The directory from which to search for
        fastq/fasta files to split.
    :param output_dir: Output directory for fastq.
    :param type: The type of sample, either Native or PCR
    :param n_bases_to_mask_tail: Number of bases to mask from the
        tail adapter (number of bases at the end of read)
    :param n_bases_to_mask_head: Number of bases to mask from the
        head adapter (number of bases at the start of read)
    :param degenerate_bases: count of N's between tail and
        head adapter (defaults to n_bases_to_mask_tail + n_bases_to_mask_head)
    :param debug_output: Output additional files helpful for debugging.
    :param edit_threshold: The threshold at which to split reads. Reads
        with edit distance below this value will be split
    :param n_replacement: Optional sequence to use as replacement
        of the masked bases. Can be used if the call is
    :param pattern: The pattern to use for matching fastq/fasta files.
    :param print_alignment: Whether to pretty-print the alignment
        at each match.
    :param print_threshold_delta: The threshold to print the
        alignment at, relative to edit_threshold
    :param threads: number of worker threads.
    :param allow_multiple_splits: Allow reads to be split more than once
    :param trim_start: How many bases to trim (mask) from the
                       beginning of the strand
    :param trim_end: How many bases to trim (mask) from the end of the strand
    """
    logger = duplex_tools.get_named_logger("SplitOnAdapters")
    logger.info(f'Duplex tools version: {duplex_tools.__version__}')
    fastxs = natsorted(list(Path(fastq_dir).rglob(pattern)), key=str)
    if output_dir is not None:
        output = Path(output_dir)
        try:
            output.mkdir()
        except FileExistsError:
            print("The output directory should not pre-exist.")
            sys.exit(1)

    targets = build_targets(
        n_bases_to_mask_head=n_bases_to_mask_head,
        n_bases_to_mask_tail=n_bases_to_mask_tail,
        degenerate_bases=degenerate_bases,
        n_replacement=n_replacement)[type]
    if edit_threshold is None:
        edit_threshold = EDIT_THRESHOLDS[type]
    edited_reads = set()
    unedited_reads = set()
    split_multiple_times = set()
    worker = functools.partial(
        process_file,
        targets=targets, output_dir=output_dir,
        debug_output=debug_output,
        edit_threshold=edit_threshold,
        print_alignment=print_alignment,
        print_threshold_delta=print_threshold_delta,
        allow_multiple_splits=allow_multiple_splits,
        trim_start=trim_start,
        trim_end=trim_end,
    )

    with ProcessPoolExecutor(max_workers=threads) as executor:
        results = executor.map(worker, fastxs)
        total_written = 0
        for edited, unedited, multi, nwritten in results:
            edited_reads.update(edited)
            unedited_reads.update(unedited)
            split_multiple_times.update(multi)
            total_written += nwritten

    with open(Path(output, 'edited.pkl'), 'wb') as handle:
        pickle.dump(edited_reads, handle)
    with open(Path(output, 'unedited.pkl'), 'wb') as handle:
        pickle.dump(unedited_reads, handle)
    if not allow_multiple_splits:
        with open(Path(output, 'split_multiple_times.pkl'), 'wb') as handle:
            pickle.dump(split_multiple_times, handle)
    nedited_reads = len(edited_reads)
    nunedited_reads = len(unedited_reads)
    n_multisplit = len(split_multiple_times)
    logger.info(f'Split {nedited_reads} reads\n'
                f'Kept {nunedited_reads} reads')
    logger.info(f'Wrote a total of {total_written} reads')
    if not allow_multiple_splits:
        logger.info(f'{n_multisplit} reads contained multiple'
                    f' adapters but we re written out as single reads '
                    f'(to split these, set --allow-multiple-splits')


def argparser():
    """Create argument parser."""
    parser = ArgumentParser(
        "Split basecalls based on adapter sequences.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        parents=[duplex_tools._log_level()], add_help=False)

    default = ' (default: %(default)s)'
    parser.add_argument(
        "fastq_dir",
        help="The directory to search for fastq/fasta files to split.")
    parser.add_argument(
        "output_dir",
        help="Output directory for fastq.")
    parser.add_argument(
        "sample_type", choices=["Native", "PCR"],
        help="Sample type.")
    parser.add_argument(
        "--pattern", default="*.fastq*",
        help="Pattern used for matching fastq/fasta files." + default)
    parser.add_argument(
        "--threads", default=None, type=int,
        help=(
            "Number of worker threads. "
            "Equal to number of logical CPUs by default."))
    parser.add_argument(
        "--n_bases_to_mask_tail", default=mask_size_default_tail, type=int,
        help=(
            "Number of bases to mask from the tail adapter "
            "(number of bases at the end of read)." + default))
    parser.add_argument(
        "--n_bases_to_mask_head", default=mask_size_default_head, type=int,
        help=(
            "Number of bases to mask from the head adapter "
            "(number of bases at the start of read)." + default))
    parser.add_argument(
        "--degenerate_bases", default=mask_size_default_N, type=int,
        help=(
            "Count of N's between tail and head adapter "
            "(defaults to n_bases_to_mask_tail + n_bases_to_mask_head)."))
    parser.add_argument(
        "--debug_output", action="store_true",
        help=(
            "Output an additional debug file to show the adapter region"))
    parser.add_argument(
        "--edit_threshold", default=None, type=int,
        help=(
            "The threshold at which to split reads. Reads with edit distance "
            "below this value will be split."))
    parser.add_argument(
        "--n_replacement", default=None, type=str,
        help="Optional sequence to use as replacement of the masked bases.")
    parser.add_argument(
        "--print_alignment", action="store_true",
        help="Pretty-print the alignment at each match.")
    parser.add_argument(
        "--print_threshold_delta", default=0, type=int,
        help="Threshold to print the alignment, "
             "relative to edit_threshold." + default)
    parser.add_argument(
        "--allow_multiple_splits", action="store_true",
        help="Whether to split reads into more than 2 new reads if there"
             " are multiple hits")
    parser.add_argument(
        "--trim_start", default=200, type=int,
        help="How many bases to trim (mask) "
             "from the beginning of the strand" + default)
    parser.add_argument(
        "--trim_end", default=200, type=int,
        help="How many bases to trim (mask) "
             "from the end of the strand" + default)
    return parser


def main(args):
    """Entry point."""
    split(
        args.fastq_dir,
        args.output_dir,
        args.sample_type,
        args.n_bases_to_mask_tail,
        args.n_bases_to_mask_head,
        args.degenerate_bases,
        args.debug_output,
        args.edit_threshold,
        args.n_replacement,
        args.pattern,
        args.print_alignment,
        args.print_threshold_delta,
        args.threads,
        args.allow_multiple_splits,
        args.trim_start,
        args.trim_end,
        )
