import tempfile
from pathlib import Path

import os

import pandas as pd

from duplex_tools.pairs_from_summary import find_pairs

import pkg_resources

seqsummary = pkg_resources.resource_filename('tests.data.summaries_for_pairing',
                                             'seqsummary.txt')
bamsummary = pkg_resources.resource_filename('tests.data.summaries_for_pairing',
                                             'dorado_unmapped.bam')


def test_pairs_from_summary():
    # Given
    outdir = 'pairs_from_seqsum'
    expected_filename = 'pair_ids.txt'
    expected_file = Path(outdir) / expected_filename
    
    # When
    find_pairs(seqsummary, outdir=outdir)
    
    # Then
    assert os.path.exists(outdir)
    assert expected_file.is_file()


def test_pairs_from_bam():
    # Given
    outdir = 'pairs_from_bam'
    expected_filename = 'pair_ids.txt'
    expected_file = Path(outdir) / expected_filename
    
    # When
    find_pairs(bamsummary, outdir=outdir)
    
    # Then
    assert os.path.exists(outdir)
    assert expected_file.is_file()


def test_pairs_are_same():
    # Given two methods
    find_pairs(bamsummary, outdir='pairs_from_bam')
    find_pairs(seqsummary, outdir='pairs_from_seqsum')
    
    
    # When comparing the pairs
    seqpairs = pd.read_csv('pairs_from_seqsum/pair_ids.txt', sep=' ',
                           names=['t', 'c'])
    bampairs = pd.read_csv('pairs_from_bam/pair_ids.txt', sep=' ',
                           names=['t', 'c'])
    ids_seq = set(seqpairs['t'] + ';' + seqpairs['c'])
    ids_bam = set(bampairs['t'] + ';' + bampairs['c'])
    missing_from_bam = ids_seq.difference(ids_bam)
    missing_from_seq = ids_bam.difference(ids_seq)
    union = ids_bam.union(ids_seq)
    intersection = ids_bam.intersection(ids_seq)
    
    # The IoU is high (> 0.9). Expect to miss some reads since
    # parent_read_id != read_id
    IoU = len(intersection) / len(union)
    print(f'Reads missing from seq: {missing_from_seq}')
    print(f'Reads missing from bam: {missing_from_bam}')
    assert IoU > 0.9
