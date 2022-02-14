#!/usr/bin/python3
"""Get answer of Life, the Universe and Everything."""

from MafBlock import MafStream
import sys
import os
import json

# Minimal size and length that a maf block must have to be processed by RNAcode.
# Absolute lower boundaries
MIN_LENGTH = 12
MIN_SIZE = 3

# Parameters for concatenation of blocks even with impeding species. Species
# will be deleted.
# Length a block must exceed so that no species will be deleted.
MIN_LENGTH_DEL = 60
# Maximal number of species that can be deleted.
MAX_DEL_SPECIES = 1

# Parameter for splitting
MAX_LEN_NO_SPLIT = 3000

# Big block size
BB_SIZE = 1000


def main():
    """Get answer of Life, the Universe and Everything."""
    maf_file_path = sys.argv[1]
    if not os.path.isfile(maf_file_path):
        print(f"Maf file {maf_file_path} does not exist!")
        sys.exit(1)

    chromosome_dir = os.path.dirname(maf_file_path)
    chromosome_name = os.path.basename(maf_file_path).split(".")[0]

    split_dir = os.path.join(os.path.split(maf_file_path)[0], "big_blocks")
    if not os.path.isdir(split_dir):
        os.mkdir(split_dir)
    else:
        for block_file in os.listdir(split_dir):
            if not os.path.isfile(os.path.join(split_dir, block_file)):
                continue
            os.remove(os.path.join(split_dir, block_file))

    maf_stream = MafStream(
        path=maf_file_path,
        min_length_del=MIN_LENGTH_DEL,
        max_del_species=MAX_DEL_SPECIES,
        min_size=MIN_SIZE,
        min_length=MIN_LENGTH,
        max_len_no_split=MAX_LEN_NO_SPLIT,
    )

    bb_num = 1
    maf_counter = 0
    block_dic = {}
    f_handle = open(os.path.join(split_dir, f"big_block_{bb_num}.maf"), "w", encoding="UTF-8")
    for maf in maf_stream.discard_stream():
        maf_counter += 1
        if maf_counter > BB_SIZE:
            bb_num += 1
            maf_counter = 0
            f_handle.close()
            f_handle = open(os.path.join(split_dir, f"big_block_{bb_num}.maf"), "w", encoding="UTF-8")
        small_target_name = f"{chromosome_name}_{maf_counter}_{bb_num}"
        block_dic[small_target_name] = maf.block_index_list
        maf.block[1][1] = small_target_name

        f_handle.write(str(maf))

    f_handle.close()
    with open(os.path.join(chromosome_dir, "block_dic.json"), "w", encoding="UTF-8") as f_handle:
        json.dump(block_dic, f_handle)


if __name__ == "__main__":
    main()


"""
    for _block_index, block in maf_stream.iterate_from(0):
        block.generate_html("/homes/biertruck/john/public_html/mview/", name="orginal_maf")
    maf_file_unzip_path = "/scr/k61san2/john/rnacode_human_CS/multiz100way/chr10_GL383545v1_alt/chr10_GL383545v1_alt.maf"
    with open(maf_file_unzip_path, "w", encoding="UTF-8") as f_handle:
    with gzip.open(maf_file_path, "rb") as maf_handle:
        f_handle.write(maf_handle.read().decode("UTF8"))

    with open("/scr/k61san2/john/rnacode_human_CS/multiz100way/chr10_GL383545v1_alt/test_new.maf", "w", encoding="UTF-8") as f_handle:
        for maf in maf_list:
            maf.sort()
            # maf.generate_html("/homes/biertruck/john/public_html/mview/", name="concat_maf_old")
            f_handle.write(str(maf))

    end_block_index, maf_list = maf_stream.concat_blocks_verbose(0, only_block=False, split=True)
    with open("/scr/k61san2/john/rnacode_human_CS/multiz100way/chr10_GL383545v1_alt/test_verbose.maf", "w", encoding="UTF-8") as f_handle:
        for maf in maf_list:
            maf.sort()
            f_handle.write(str(maf))
    end_block_index, maf_list = maf_stream.concat_blocks(0, only_block=False, split=True)
    with open("/scr/k61san2/john/rnacode_human_CS/multiz100way/chr10_GL383545v1_alt/test.maf", "w", encoding="UTF-8") as f_handle:
        for maf in maf_list:
            maf.sort()
            f_handle.write(str(maf))
"""
