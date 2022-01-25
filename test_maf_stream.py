#!/usr/bin/python3
"""Get answer of Life, the Universe and Everything."""

from MafBlock import MafStream
from MafBlock import MafBlock

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


def concat_split(maf_file_path):
    """Get answer of Life, the Universe and Everything."""
    maf = MafBlock()
    maf_list = []

    maf_stream = MafStream(maf_file_path)
    for block_index, next_maf in maf_stream:
        next_maf.add_suffix(f"block-index_{block_index}")
        # Only first iteration
        if maf.is_empty():
            maf = next_maf
            continue
        ret_val = maf.concat(next_maf)
        # Return value is zero if concatenation was successful
        if ret_val == 0:
            pass
        # The returned negative int says how many species imped a concatenation.
        elif ret_val >= -MAX_DEL_SPECIES:
            # If the current block or the next block are below threshold force
            # concatenation even with deleting species.
            try:
                if len(maf) < MIN_LENGTH_DEL or len(maf_stream.concat_blocks(block_index + 1)) < MIN_LENGTH_DEL:
                    maf.concat(next_maf, max_del=MAX_DEL_SPECIES)
            except IndexError as e:
                # Catch EOF
                if str(e) == "block index out of range":
                    ret_val = -MAX_DEL_SPECIES - 1
                else:
                    raise
            # If both blocks are sufficiently big keep current block and proceed
            # with concatenation with next block.
            else:
                ret_val = -MAX_DEL_SPECIES - 1

        # Here concatenation ends. Either to many impeding species or impeding
        # species is below MAX_DEL_SPECIES but both maf blocks, the current and
        # the next are big enough.
        if ret_val < -MAX_DEL_SPECIES:
            # Catches lower bound size
            if (
                maf.size() - 1 < MIN_SIZE
                and len(maf.block[1][-1].replace("-", "")) < MIN_LENGTH
            ):
                maf = next_maf
                continue
            # Cut blocks
            maf_list += maf.preprocess_block(MAX_LEN_NO_SPLIT)
            maf = next_maf


def main():
    """Get answer of Life, the Universe and Everything."""
    maf_file_path = "/scr/k61san2/john/rnacode_human_CS/multiz100way/chr10_GL383545v1_alt/chr10_GL383545v1_alt.maf.gz"
    maf_stream = MafStream(
        path=maf_file_path,
        min_length_del=MIN_LENGTH_DEL,
        max_del_species=MAX_DEL_SPECIES,
        min_size=MIN_SIZE,
        min_length=MIN_LENGTH,
        max_len_no_split=MAX_LEN_NO_SPLIT,
    )
    maf_list = maf_stream.concat_blocks(0, only_block=False, split=True)
    print(len(maf_list))


if __name__ == "__main__":
    main()
