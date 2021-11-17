import io
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

import os
import pkg_resources
import shutil
from duplex_tools.split_on_adapter import split

def test_split_on_adapter():
    # Given some data
    fastq_dir = pkg_resources.resource_filename('tests.data', 'fastq_200-th-200')
    filename_in = pkg_resources.resource_filename('tests.data.fastq_200-th-200', '200bases-tailhead-200bases.fastq')
    expected_file = Path(filename_in.replace('.fastq', '_split.fastq.gz')).name

    # When running the script on it
    output = io.StringIO()

    # Make sure the dir does not exist before, we want to check it's created
    dir = 'dir_with_split_reads'
    shutil.rmtree(dir, ignore_errors=True)
    with redirect_stdout(output):
        split(fastq_dir, output_dir=dir, pattern='*.fastq')

    # Then (1) there are files created in the expected locations
    assert os.path.exists(dir)
    assert (Path(dir) / expected_file).is_file()

    # Then (2) stdout from script contains the number of split reads
    print(output.getvalue())
    assert "Split 1 reads" in output.getvalue()

def test_split_gzipped_and_nonzipped_by_default():
    # Given some data
    fastq_dir = pkg_resources.resource_filename('tests.data',
                                                'fastq_200-th-200')
    output = io.StringIO()
    tmpdir = next(tempfile._get_candidate_names())

    with redirect_stdout(output):
        split(fastq_dir, output_dir=tmpdir)

    shutil.rmtree(tmpdir)

    assert 'Split 2 reads' in output.getvalue()
