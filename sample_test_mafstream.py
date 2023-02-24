#!/usr/bin/python3
"""Auxilary script to get some maf blocks by range or length."""

from MafBlock import MafStream


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
    maf_file_path = "/scr/k61san2/john/rnacode_human_CS/multiz100way/chr1/chr1.maf.gz"
    # maf_file_path = "/scr/k61san2/john/rnacode_human_CS/multiz100way/chr1/big_blocks/big_block_1.maf"
    maf_file_path = "./test.maf"

    maf_stream = MafStream(
        path=maf_file_path,
        min_length_del=MIN_LENGTH_DEL,
        max_del_species=MAX_DEL_SPECIES,
        min_size=MIN_SIZE,
        min_length=MIN_LENGTH,
        max_len_no_split=MAX_LEN_NO_SPLIT,
    )

    # index_range=(6954439, 6954463)
    # for maf in maf_stream.iterate_from(index_range[0]):
    #     if max(maf.block_index_list) > index_range[1]:
    #         print(maf)
    #         break
    #     print(maf)
    # exit()
    for maf in maf_stream.discard_stream():
        print(maf)
    exit()
    # for maf in maf_stream:
    #     print(maf)
    #     input("Wait")

    for maf in maf_stream:
        pass
    print(maf)

    # for maf in maf_stream.discard_stream():
    #     print(maf)
    #     print(maf.size())

    #     print(maf.size() - 1 < maf_stream.min_size)
    #     print(maf.len_no_gaps() < maf_stream.min_length)
    #     input("Wait")


if __name__ == "__main__":
    main()
