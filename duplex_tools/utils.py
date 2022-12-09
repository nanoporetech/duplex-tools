"""Utilities for duplex-tools."""
import pysam


def contains_references(pysam_bam: pysam.AlignmentFile):
    """Check whether bam contains any references."""
    return pysam_bam.header.nreferences > 0


def is_ubam(pysam_bam: pysam.AlignmentFile):
    """Check whether a bam is uBAM."""
    return not contains_references(pysam_bam)
