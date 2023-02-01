"""Methods for splitting reads using a sam/bam with moves and a pod5.

1. One method for finding split points.
2. One method for splitting raw data into new reads based on 1.
"""
from collections import defaultdict
from concurrent.futures import as_completed, ProcessPoolExecutor
from functools import partial
from pathlib import Path
import random
import uuid

from matplotlib import pyplot as plt
from more_itertools import chunked
import pandas as pd
import pod5 as p5
from pod5 import EndReason, EndReasonEnum, Read, Reader
import pysam
from tqdm import tqdm

import duplex_tools
from duplex_tools.split_pairs_utils import split


def get_split_points(
        input_dorado_xam,
        threads=-1,
        max_reads: int = None,
        chunk_size=5000,
        match_threshold=0.8,
        left_midpoint_threshold=0.45,
        right_midpoint_threshold=0.55
) -> dict:
    """Detect locations where reads in the Dorado SAM/BAM file should be split.

    :param left_midpoint_threshold: Require the midpoint to be -> of here.
    :param right_midpoint_threshold: Require the midpoint to be <- of here.
    :param match_threshold: Require at least this fraction of the template
           to match the complement
    :param input_dorado_xam: The path to the input Dorado SAM/BAM file.
    :param threads: The number of threads to use.
    :param max_reads: The maximum number of reads to process.
    :param chunk_size: The number of reads to process in one batch
                       (use lower value for more frequent logging)
    :return: A dictionary containing the split locations for each read.
    """
    logger = duplex_tools.get_named_logger("SplitPairs")

    splitter = partial(split,
                       match_threshold=match_threshold,
                       left_midpoint_threshold=left_midpoint_threshold,
                       right_midpoint_threshold=right_midpoint_threshold
                       )
    split_locations = {}
    finished = False
    counter = defaultdict(int)
    with ProcessPoolExecutor(threads) as pool:
        with pysam.AlignmentFile(input_dorado_xam, check_sq=False) as f:
            it = f.fetch(until_eof=True)
            iterator = chunked(it, chunk_size)

            for idx, chunk in enumerate(iterator):
                if finished:
                    break
                futures = []

                for read in chunk:
                    seq = read.query_sequence
                    elem = (
                        read.qname,
                        seq,
                        read.get_tag("mv"),
                        read.get_tag("ts"),
                        read.get_tag("ns"),  # Sample count (For signal end)
                    )
                    future = pool.submit(splitter, elem)
                    futures.append(future)

                counter['assessed'] += len(futures)
                for f in as_completed(futures):
                    if f.result() is not None:
                        split_locations.update(f.result())

                counter['tosplit'] = len(split_locations.keys())

                if max_reads is not None and len(split_locations.keys()) > \
                        max_reads:
                    finished = True
                assessed = counter['assessed']
                nsplit = counter['tosplit']
                logger.info(
                    "Split/Processed reads:"
                    f"{nsplit:.0f}/{assessed:.0f}"
                    f" ({100 * nsplit / assessed:.2f}%)")
    logger.info("Finished finding breakpoints.")
    return split_locations


def split_pod5(
        input_pod5_dir, split_locations, new_pod5_dir, force_overwrite=False,
        debug_dir=None
):
    """Split a pod5 given pre-determined locations to split them.

    :param debug_dir: A directory in which to write plots of raw signal
    :param input_pod5_dir: The directory containing the pod5 file(s) to be
                           split.
    :param split_locations: A dictionary of read IDs and their corresponding
                            split locations.
    :param new_pod5_dir: The directory where the split pod5 file will be saved.
    :param force_overwrite: If True, any existing files in `new_pod5_dir` will
                            be overwritten.
    :return: None

    """
    logger = duplex_tools.get_named_logger("SplitPairs")

    Path(new_pod5_dir).mkdir(exist_ok=True, parents=True)

    ids = list(split_locations.keys())
    if debug_dir is not None:
        Path(debug_dir).mkdir(exist_ok=True, parents=True)

    pod5_dir_input = input_pod5_dir
    for pod5 in Path(pod5_dir_input).rglob("*.pod5"):
        read_ids = []
        logger.info(f"Splitting {pod5}")
        output_stem = pod5.stem + "_split_duplex"
        p5_file = f"{str(new_pod5_dir)}/{output_stem}.pod5"
        if force_overwrite:
            Path(p5_file).unlink(missing_ok=True)
        with p5.Writer(p5_file) as writer:
            for read in tqdm(Reader(pod5).reads(selection=ids,
                                                missing_ok=True)):
                read_obj = read.to_read()
                data = vars(read_obj)
                rd = random.Random()
                rd.seed(read.read_id)
                id_left = uuid.UUID(int=rd.getrandbits(128), version=4)
                id_right = uuid.UUID(int=rd.getrandbits(128), version=4)
                read_id_map = {
                    "read_id": read.read_id,
                    "read_id_left": id_left,
                    "read_id_right": id_right,
                }
                read_ids.append(read_id_map)

                at = split_locations[str(read.read_id)]

                # Grab the data for the template read
                signal_orig = read.signal
                signal_left = data["signal"][at["left"][0]: at["left"][1]]
                signal_right = data["signal"][at["right"][0]: at["right"][1]]
                data["signal"] = signal_left
                data["read_id"] = read_id_map["read_id_left"]
                data["end_reason"] = EndReason(EndReasonEnum.UNKNOWN, False)
                left = Read(**data)

                # Grab the data for the complement read
                data["signal"] = signal_right
                data["read_id"] = read_id_map["read_id_right"]
                data["end_reason"] = EndReason(EndReasonEnum.UNKNOWN, False)
                right = Read(**data)

                # Write the read objects
                writer.add_read(left)
                writer.add_read(right)
                if debug_dir is not None:
                    _, axes = plt.subplots(3, 1, figsize=(23, 6))
                    for idx, (data, ax) in enumerate(zip([signal_orig,
                                                          signal_left,
                                                          signal_right],
                                                         axes)):
                        ax.plot(data, linewidth=0.05)
                        if idx == 0:
                            ax.axvline(at["left"][1], linewidth=0.2,
                                       color='red')
                    plt.savefig(f'{debug_dir}/{read.read_id}.png')
                    plt.close()

            df = pd.DataFrame(read_ids)
            df.to_csv(
                f"{new_pod5_dir}/{output_stem}.txt", sep="\t", index=False
            )
            try:
                df[["read_id_left", "read_id_right"]].to_csv(
                    f"{new_pod5_dir}/{output_stem}_pair_ids.txt",
                    sep=" ",
                    index=False,
                )
                logger.info(f"Created {len(read_ids)} new pairs")
            except KeyError:
                logger.info("No pairs created")
