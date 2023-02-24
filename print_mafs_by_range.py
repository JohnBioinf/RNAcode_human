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

maf_stream = MafStream(
    path="/scr/k61san2/john/rnacode_human_CS/multiz100way/chr14/chr14.maf.gz",
    min_length_del=MIN_LENGTH_DEL,
    max_del_species=MAX_DEL_SPECIES,
    min_size=MIN_SIZE,
    min_length=MIN_LENGTH,
    max_len_no_split=MAX_LEN_NO_SPLIT,
)

maf_stream = MafStream(
    path="./test.maf",
    min_length_del=MIN_LENGTH_DEL,
    max_del_species=MAX_DEL_SPECIES,
    min_size=MIN_SIZE,
    min_length=MIN_LENGTH,
    max_len_no_split=MAX_LEN_NO_SPLIT,
)

start = 1707553
end = 1707558
index_range = (1707553, 1707560)
position = (16022637, 16022669)
position = (64159073, 64159073 + 50)

print("Concat from index range")
for maf in maf_stream.split_stream():
    print("Start")
    print(maf)
    print(maf.coordinates())
exit()

print("Print blocks by index")
for maf in maf_stream.iterate_from(0):
    print(maf)
    if maf.block_index_list[0] > index_range[1]:
        break
exit()

print("Concat from index range")
for maf in maf_stream.split_stream(index_range=index_range):
    print("Start")
    print(maf)
    print(maf.coordinates())
exit()

print("Print blocks by index")
for maf in maf_stream.iterate_from(index_range[0]):
    print(maf)
    if maf.block_index_list[0] > index_range[1]:
        break
exit()

print("Concat from position")
for maf in maf_stream.discard_stream(position=position):
    print("Start")
    print(maf)
    print(maf.coordinates())

exit()

exit()
