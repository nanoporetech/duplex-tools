import io
from contextlib import redirect_stdout
from pathlib import Path

import os
import pkg_resources
import shutil
from duplex_tools.split_on_adapter import split

def test_read_fillet():
    # Given some data
    fastq_dir = pkg_resources.resource_filename('tests.data', 'fastq_200-th-200')
    filename_in = pkg_resources.resource_filename('tests.data.fastq_200-th-200', '200bases-tailhead-200bases.fastq')
    expected_file = Path(filename_in.replace('.fastq','_split.fastq.gz')).name

    # When running the script on it
    output = io.StringIO()
    dir = 'dir_with_split_reads'
    shutil.rmtree(dir, ignore_errors=True)
    with redirect_stdout(output):
        split(fastq_dir, output_dir=dir, pattern='*.fastq')

    # Then (1) there are files created in the expected locations
    assert os.path.exists(dir)
    assert (Path(dir) / expected_file).is_file()

    # Then (2) the stdout from the script has information about the number of split reads
    print(output.getvalue())
    assert "Split 1 reads" in output.getvalue()
