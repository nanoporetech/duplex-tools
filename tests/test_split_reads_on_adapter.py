import io
import tempfile
from contextlib import redirect_stderr
from pathlib import Path

import os

import duplex_tools
import pkg_resources
import shutil
import logging
from duplex_tools.split_on_adapter import split


def test_split_on_adapter(caplog):
    caplog.set_level(logging.INFO)
    # Given some data
    fastq_dir = pkg_resources.resource_filename('tests.data', 'fastq_200-th-200')
    filename_in = pkg_resources.resource_filename('tests.data.fastq_200-th-200', '200bases-tailhead-200bases.fastq')
    expected_file = Path(filename_in.replace('.fastq', '_split.fastq.gz')).name

    # When running the script on it
    output = io.StringIO()

    # Make sure the dir does not exist before, we want to check it's created
    dir = 'dir_with_split_reads'
    shutil.rmtree(dir, ignore_errors=True)

    split(fastq_dir, output_dir=dir, pattern='*.fastq')

    # Then (1) there are files created in the expected locations
    assert os.path.exists(dir)
    assert (Path(dir) / expected_file).is_file()

    # Then (2) stdout from script contains the number of split reads
    print(output.getvalue())
    assert "Split 1 reads" in ''.join(caplog.text)


def test_split_gzipped_and_nonzipped_by_default(caplog):
    caplog.set_level(logging.INFO)
    # Given some data
    fastq_dir = pkg_resources.resource_filename('tests.data',
                                                'fastq_200-th-200')
    tmpdir = next(tempfile._get_candidate_names())

    split(fastq_dir, output_dir=tmpdir)

    shutil.rmtree(tmpdir)

    assert 'Split 2 reads' in ''.join(caplog.text)


def test_split_multiple_on_adapter(caplog):
    caplog.set_level(logging.INFO)
    # Given some data
    fastq_dir = pkg_resources.resource_filename('tests.data', 'fastq_200-th-200-th-200')
    filename_in = pkg_resources.resource_filename('tests.data.fastq_200-th-200-th-200',
                                                  '200bases-tailhead-200bases-tailhead-200bases.fastq')
    expected_file = Path(filename_in.replace('.fastq', '_split.fastq.gz')).name

    # When running the script on it

    logger = duplex_tools.get_named_logger("SplitOnAdapters")
    logger.info('Test')
    # Make sure the dir does not exist before, we want to check it's created
    dir = 'dir_with_split_reads'
    shutil.rmtree(dir, ignore_errors=True)
    split(fastq_dir, output_dir=dir, pattern='*.fastq', allow_multiple_splits=True, trim_end=20)

    # Then (1) there are files created in the expected locations
    assert os.path.exists(dir)
    assert (Path(dir) / expected_file).is_file()

    # Then (2) stdout from script contains the number of split reads
    assert "Split 1 reads" in ''.join(caplog.text)
