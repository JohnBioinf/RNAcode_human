#!/usr/bin/python3
"""Remove all erronous maf blocks.

In the computation and anlysis of the compelete concatenatd WGA a bug crippled
some maf blocks with to many gaps. These blocks are of low quality and thus can
be simply discarded with only little statiscial impact."""


import json

with open("./parameters_local.json", "r", encoding="UTF-8") as FILE_HANDLE:
    PARAMETERS = json.load(FILE_HANDLE)

DATA_DIR = PARAMETERS["DATA_DIR"]
DATA_DIR_OLD = PARAMETERS["DATA_DIR_OLD"]
MULTIZ100WAY_DIR = DATA_DIR + "/multiz100way/"
MULTIZ100WAY_DIR_OLD = DATA_DIR_OLD + "/multiz100way/"

with open("./fuck_up.lst", "r", encoding="UTF-8") as f_handle:
    FUCK_UP_LIST = f_handle.read().split("\n")


def clean(genome_alignment_dir, remove_from_list=True):
    """Get answer of Life, the Universe and Everything."""
    org_bed_annotation = genome_alignment_dir + "/RNAcode.bed"
    clean_bed_annotation = genome_alignment_dir + "/RNAcode_clean.bed"

    with open(org_bed_annotation, "r", encoding="UTF-8") as f_handle:
        anno = [line.split() for line in f_handle.read().split("\n")]

    with open(clean_bed_annotation, "w", encoding="UTF-8") as f_handle:
        for line in anno:
            if line == []:
                continue

            if "_" in line[0]:
                continue

            if int(line[2]) - int(line[1]) < 10:
                print(line)
                continue

            if remove_from_list:
                is_a_fuck_up = False
                for block in line[3].split(","):
                    hss_id = block.split("-")[1]
                    if hss_id in FUCK_UP_LIST:
                        is_a_fuck_up = True
                        break
                if is_a_fuck_up:
                    continue

            f_handle.write(" ".join(line) + "\n")
        f_handle.write("\n")


def main():
    """Get answer of Life, the Universe and Everything."""
    clean(MULTIZ100WAY_DIR, remove_from_list=True)
    # clean(MULTIZ100WAY_DIR_OLD, remove_from_list=False)


if __name__ == "__main__":
    main()
