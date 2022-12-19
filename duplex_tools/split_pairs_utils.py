"""Utilities for split_pairs_steps.py."""
import mappy as mp
import numpy as np


def map_to_itself(seq):
    """Map a sequence to its reverse complement."""
    aligner = mp.Aligner(seq=seq, extra_flags=0x200000)
    for match in aligner.map(seq):
        yield match


def detect_duplex_location(seq) -> tuple:
    """
    Find an optional duplex self-alignment.

    For a duplex self-alignment to be valid,
    the ref and query starts and ends have to be identical

    The midpoint of such an alignment is the point at which we want to split
    the read
    """
    for ali in map_to_itself(seq):
        if (ali.r_st == ali.q_st) & (ali.r_en == ali.q_en):

            midpoint = int((ali.r_st + ali.r_en) / 2)
            if midpoint < len(seq) * 0.1:
                continue
            if midpoint > len(seq) * 0.9:
                continue
            tempcomp_lenratio = (len(seq) - midpoint) / (midpoint - ali.r_st)
            if tempcomp_lenratio < 3:
                return (0, ali.r_st, midpoint, ali.r_en, len(seq))


def split(read,
          match_threshold=0.8,
          left_midpoint_threshold=0.45,
          right_midpoint_threshold=0.55
          ) -> dict:
    """Map reads to themselves.

    Emit the sequences as new reads.


    >>> split(['TTTTTTTTTTTTGGAGATCTCCAAAAAAA'])
                                ^^

    ('TTTTTTTTTTTTGGAGA', (0, 17))
    ('TCTCCAAAAAAA', (17, 29))

    :param match_threshold: The fraction of the strand required to self-align
    :param left_midpoint_threshold: Midpoint must be after this fraction
                                    of bases
    :param right_midpoint_threshold: Midpoint must be before this fraction
                                     of bases
    :param read: A tuple with (read_id, sequence, moves, trimmed_samples)

    """
    read_id, sequence, mv, ts, nsamples = read
    maps_to_self, coordinates = selfmap_mini(sequence)
    # First element in the move array contains the stride
    stride = mv[0]
    # Rest of elements contain the actual moves
    mv = mv[1:]
    if not maps_to_self:
        return None
    mv_cs = np.cumsum(mv)

    left_alistart = coordinates["1_alistart"]
    left_end = coordinates["2_midpoint"]
    right_start = coordinates["2_midpoint"]

    match_frac = (left_end - left_alistart)/(left_end)
    if match_frac < match_threshold:
        return None

    if not (left_midpoint_threshold*len(sequence) <
            left_end <
            right_midpoint_threshold*len(sequence)):
        return None

    left_start_signal = 0
    left_end_signal = np.searchsorted(mv_cs, left_end) * stride + ts
    right_start_signal = np.searchsorted(mv_cs, right_start) * stride + ts
    right_end_signal = nsamples

    elem = {
        read_id: {
            "left": (left_start_signal, left_end_signal),
            "right": (right_start_signal, right_end_signal),
        }
    }
    return elem


def selfmap_mini(seq):
    """Map a sequence onto itself and get split points."""
    duploc = detect_duplex_location(seq)
    if duploc is None:
        return False, {}
    _0refstart, _1alistart, _2midpoint, _3aliend, _4seqend = duploc
    coordinates = {
        "0_refstart": _0refstart,
        "1_alistart": _1alistart,
        "2_midpoint": _2midpoint,
        "3_aliend": _3aliend,
        "4_seqend": _4seqend,
    }
    return True, coordinates
