#!/usr/bin/python3
"""Calculates the FDR for the recall."""

import numpy as np
import json
import os
from bisect import bisect_right
import sys


def sum_scores(annotation_file_path, hss_score_dic, out_file_path):
    """Get answer of Life, the Universe and Everything."""
    p_val_list = []

    with open(annotation_file_path, "r", encoding="UTF-8") as f_handle:
        for line in f_handle:
            # legacy
            # hss_id = f"{line.split()[0]}_{line.split()[3]}"
            hss_id = line.split()[3]
            p_val_list.append(hss_score_dic[hss_id])

    # uniq_p_val = list(set(p_val_list))
    # uniq_p_val.sort()
    p_val_list.sort()

    with open(out_file_path, "w", encoding="UTF-8") as f_handle:
        # for i in reversed(np.arange(1, 17, 0.000001)):
        for p_val in np.logspace(-1, -17, num=100):
            # p_val = pow(10, -i)
            f_handle.write(f"{p_val}\t{bisect_right(p_val_list, p_val)}\n")


def calc_fdr_recall(genome_alignment_dir):
    """Get answer of Life, the Universe and Everything."""
    hss_score_dic = {}
    for chromosome in os.listdir(genome_alignment_dir):
        if not os.path.isdir(genome_alignment_dir + "/" + chromosome):
            continue
        chromosome_dir_path = f"{genome_alignment_dir}/{chromosome}/"
        hss_score_dic_path = f"{chromosome_dir_path}/hss_score_dic.json"
        with open(hss_score_dic_path, "r", encoding="UTF-8") as f_handle:
            hss_score_dic.update(json.load(f_handle))

    annotation_file_list = [
        "RNAcode_overlap",
        "RNAcode_no_overlap",
    ]

    for annotation_file in annotation_file_list:
        annotation_file_path = f"{genome_alignment_dir}/{annotation_file}.bed"
        result_file_path = f"{genome_alignment_dir}/{annotation_file}.tsv"
        sum_scores(annotation_file_path, hss_score_dic, result_file_path)


def main():
    """Get answer of Life, the Universe and Everything."""
    genome_alignment_dir = sys.argv[1]
    calc_fdr_recall(genome_alignment_dir)


if __name__ == "__main__":
    main()
