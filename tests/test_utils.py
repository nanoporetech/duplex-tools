import pkg_resources
import pysam

from duplex_tools.utils import is_ubam

mapped = pkg_resources.resource_filename('tests.data.summaries_for_pairing',
                                             'dorado_mapped.bam')
unmapped = pkg_resources.resource_filename('tests.data.summaries_for_pairing',
                                             'dorado_unmapped.bam')


def test_is_ubam():
    assert not is_ubam(pysam.AlignmentFile(mapped, check_sq=False))
    assert is_ubam(pysam.AlignmentFile(unmapped, check_sq=False))