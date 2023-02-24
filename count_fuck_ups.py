#!/usr/bin/python3
"""List all erronous maf blocks.

In the computation and anlysis of the compelete concatenatd WGA a bug crippled
some maf blocks with to many gaps. These blocks are of low quality and thus can
be simply discarded with only little statiscial impact."""


import json
from glob import glob
from MafBlock import MafStream
import os

with open("./parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
MULTIZ100WAY_DIR = DATA_DIR + "/multiz100way/"

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

maf_count = 0
fuck_up_count = 0
fuck_up_list = []

def get_chrom_size(maf_file_path):
    maf_stream = MafStream(
        path=maf_file_path,
        min_length_del=MIN_LENGTH_DEL,
        max_del_species=MAX_DEL_SPECIES,
        min_size=MIN_SIZE,
        min_length=MIN_LENGTH,
        max_len_no_split=MAX_LEN_NO_SPLIT,
    )
    for maf in maf_stream:
        return maf.block[1][5]


for chromosome in glob(MULTIZ100WAY_DIR + "/*"):
    if not os.path.isdir(chromosome):
        continue

    print(chromosome)
    correct_chromsize = get_chrom_size(glob(chromosome + "/*maf.gz")[0])
    
    for big_block in glob(chromosome + "/big_blocks/*maf"):
        maf_stream = MafStream(
            path=big_block,
            min_length_del=MIN_LENGTH_DEL,
            max_del_species=MAX_DEL_SPECIES,
            min_size=MIN_SIZE,
            min_length=MIN_LENGTH,
            max_len_no_split=MAX_LEN_NO_SPLIT,
        )
        for maf in maf_stream.iterate_from(0):
            maf_count += 1
            if maf.block[1][5] != correct_chromsize:
                fuck_up_count += 1
                fuck_up_list.append(maf.block[1][1])

print(fuck_up_count)
print(maf_count)

with open("./fuck_up.lst", "w", encoding="UTF-8") as f_handle:
    f_handle.write("\n".join(fuck_up_list) + "\n")

print(maf_count/fuck_up_count)
