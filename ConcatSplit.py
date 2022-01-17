#!/usr/bin/python3
"""Concatinate each maf block with the next in a multi maf file and afterwards cut them up with overlapping windows.

Core functionality can be found in MafBlock module.

:param str chromosome_dir_path: Path to chromome allignment.
:param str chromosom: Name of the chromosome.
:param int min_size: Minimal size of a maf block.
:param int min_length: Minimal length of a maf block.
:param int max_len_no_split: Maximal length of a maf block before it will be splitted.
:param int num_cpus: Number of process spwaned in parallel.
"""


import sys
import os
from glob import glob
import gzip

from MafBlock import MafBlock


def get_arguments():
    """Read arguments and set them global."""
    global CHROMOSOME_DIR_PATH, CHROMOSOME, MIN_SIZE, MIN_LENGTH
    global MAX_LEN_NO_SPLIT, NUM_CPUS
    CHROMOSOME_DIR_PATH = sys.argv[1]
    CHROMOSOME = sys.argv[2]
    MIN_SIZE = int(sys.argv[3])
    MIN_LENGTH = int(sys.argv[4])
    MAX_LEN_NO_SPLIT = int(sys.argv[5])
    NUM_CPUS = int(sys.argv[6])


def concat_and_split():
    """Build big maf-blocks."""
    chro_maf_gz_file_path = f"{CHROMOSOME_DIR_PATH}/{CHROMOSOME}.maf.gz"
    concat_maf_path = f"{CHROMOSOME_DIR_PATH}/concat_split"
    if not os.path.isdir(concat_maf_path):
        os.mkdir(concat_maf_path)
    else:
        for a_file in glob(f"{concat_maf_path}/*"):
            os.remove(a_file)

    block_index = 0
    maf = MafBlock()
    next_maf = MafBlock()
    big_block = []

    with gzip.open(chro_maf_gz_file_path, "rb") as file_handle_gz:
        for line in file_handle_gz:
            line = line.decode("UTF8")
            if line != "\n":
                next_maf.add(line)
            else:
                block_index += 1
                next_maf.add_suffix(f"_block_index-{block_index}")
                if maf.is_empty():
                    maf = next_maf
                    next_maf = MafBlock()
                    continue
                try:
                    maf.concat(next_maf)
                except ValueError as e:
                    if str(e) == "Maf block self and other are not continuous":
                        if (
                            maf.size() - 1 < MIN_SIZE
                            and len(maf.block[1][-1].replace("-", "")) < MIN_LENGTH
                        ):
                            maf = next_maf
                            next_maf = MafBlock()
                            continue

                        big_block += maf.preprocess_block(MAX_LEN_NO_SPLIT)
                maf = next_maf
                next_maf = MafBlock()

    # find appropriate bin size
    for bin_size in [1000, 100, 10, 1]:
        big_block_size = len(big_block) / bin_size
        if big_block_size > NUM_CPUS:
            break
    # print(f"Bin size is {bin_size}")

    for i, index in enumerate(range(0, len(big_block), bin_size)):
        with open(
            f"{concat_maf_path}/big_block_{i}.maf", "w", encoding="UTF-8"
        ) as f_handle:
            for maf in big_block[index : index + bin_size]:
                f_handle.write(str(maf))


def main():
    """Do it all."""
    get_arguments()
    concat_and_split()


if __name__ == "__main__":
    main()
