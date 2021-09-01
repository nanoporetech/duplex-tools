"""Assessment of read_fillet results."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import pickle

import numpy as np
import pandas as pd

import duplex_tools


def assess(
        seqkit_stats_nosecondary, edited_reads, unedited_reads,
        split_multiple_times, suffix=None):
    """Run assessment."""
    with open(edited_reads, 'rb') as handle:
        edited_reads = pickle.load(handle)
    with open(unedited_reads, 'rb') as handle:
        unedited_reads = pickle.load(handle)
    with open(split_multiple_times, 'rb') as handle:
        split_multiple_times = pickle.load(handle)

    txt = pd.read_csv(
        seqkit_stats_nosecondary,
        sep='\t')
    print(f'{seqkit_stats_nosecondary} contains {len(txt)} reads')

    expected_read_ids = (
        set(edited_reads)
        .union(unedited_reads)
        .union(split_multiple_times))
    txt = txt[txt['Read'].isin(expected_read_ids)]
    print(f'Using {len(txt)} reads for assessment')

    txt['qstart'] = np.where(
        txt['Strand'] == 1, txt['LeftClip'], txt['RightClip'])
    txt['qend'] = np.where(
        txt['Strand'] == -1,
        txt['ReadLen'] - txt['LeftClip'], txt['ReadLen'] - txt['RightClip'])

    txt['read_count'] = txt.groupby('Read')['Ref'].transform('count')
    txt['bases_between_reads'] = np.where(
        txt['read_count'] == 2,
        (txt.groupby('Read')['qstart'].transform('max')
            - txt.groupby('Read')['qend'].transform('min')),
        0)
    txt['bases_between_alignments_start'] = np.where(
        txt['read_count'] == 2,
        (txt.groupby('Read')['Pos'].transform('max')
            - txt.groupby('Read')['Pos'].transform('min')),
        int(1e8))
    txt['bases_between_alignments_end'] = np.where(
        txt['read_count'] == 2,
        (txt.groupby('Read')['Pos'].transform('max')
            - txt.groupby('Read')['Pos'].transform('min')),
        int(1e8))
    txt['bases_between_alignments_min'] = np.where(
        ((txt['read_count'] == 2)
            & (txt.groupby('Read')['Ref'].transform('nunique') == 1)),
        txt[['bases_between_alignments_start',
            'bases_between_alignments_end']].min(axis=1),
        int(1e8))

    txt['expected_class'] = None
    txt.loc[
        (txt['ReadCov'] > 95) & (txt['read_count'] == 1),
        'expected_class'] = 'single_alignment_95%cov'
    txt.loc[
        (txt['bases_between_reads'].between(20, 160)
            & (txt['read_count'] == 2)),
        'expected_class'] = 'disjoint_with_gap'
    txt.loc[
        txt['bases_between_reads'].lt(20) & (txt['read_count'] == 2),
        'expected_class'] = 'disjoint_without_gap'
    txt.loc[
        txt['bases_between_alignments_min'] < 10000,
        'expected_class'] = 'overlapping'
    txt.loc[
        txt['read_count'] > 2, 'expected_class'] = 'Read_gt-2-supplementary'

    txt = txt.query('expected_class != "overlapping"')

    txt['split_class'] = "None"
    txt.loc[txt['Read'].isin(edited_reads), 'split_class'] = 'Read split'
    txt.loc[txt['Read'].isin(unedited_reads), 'split_class'] = 'Read not split'
    txt.loc[
        txt['Read'].isin(split_multiple_times),
        'split_class'] = 'Split more than once'

    txt_only_once_twice = txt.query(
        '(split_class!= "Split more than once") & (split_class != "None") ')

    crosstabbed = pd.crosstab(
        txt_only_once_twice['split_class'],
        txt_only_once_twice['expected_class'])
    crosstabbed_fraction = crosstabbed.apply(lambda r: r / r.sum(), axis=0)
    crosstabbed.assign(label=suffix).to_csv(
        f'read_splitting_assessment_{suffix}.txt', sep='\t')
    crosstabbed_fraction.assign(label=suffix).to_csv(
        f'read_splitting_assessment_{suffix}_fraction.txt', sep='\t')

    print('Amount of data excluded')
    print(1 - len(txt_only_once_twice) / len(txt))


def argparser():
    """Create argument parser."""
    parser = ArgumentParser(
        (
            "Assess the results of the split_on_adapter command. "
            "This tool is not intended for public consumption."),
        formatter_class=ArgumentDefaultsHelpFormatter,
        parents=[duplex_tools._log_level()], add_help=False)
    parser.add_argument(
        "seqkit_stats_nosecondary")
    parser.add_argument(
        "edited_reads")
    parser.add_argument(
        "unedited_reads")
    parser.add_argument(
        "split_multiple_times")
    parser.add_argument(
        "--suffix")
    return parser


def main(args):
    """Entry point."""
    assess(
        args.seqkit_stats_nosecondary,
        args.edited_reads,
        args.unedited_reads,
        args.split_multiple_times,
        args.suffix)
