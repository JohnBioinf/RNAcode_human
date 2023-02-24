#!/usr/bin/python3
"""Calculates the recall rate and FDR for RNAcode annotation.

The calculaten works under the assumption that the current annoation is the ground truth.
"""

import json
import os
from bisect import bisect_right
import sys

import numpy as np

from IntervalFrameTree import IntervalFrameTree


with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)
DATA_DIR = parameters["DATA_DIR"]
HG38_ANNOTATION_WEB_FTP = parameters["HG38_ANNOTATION_WEB_FTP"]
ANNOTATION_FILE = HG38_ANNOTATION_WEB_FTP.split("/")[-1].replace(".gz", "")
ANNOTATION_FILE_MOD = ANNOTATION_FILE.replace(".gtf", ".cds.gtf")
MULTIZ100WAY_DIR = DATA_DIR + "/multiz100way/"
MULTIZ20WAY_DIR = DATA_DIR + "/multiz20way/"

# def recall_rate(gff_index_dic, rnacode_index_dic, genome_alignment_dir):
#     """Calculate how many protein genes are recalled."""
#     recall = {}
#     protein_id = None
#     collecting_exons = False
#     with open(f"{DATA_DIR}/{ANNOTATION_FILE}", "r", encoding="UTF-8") as f_handle:
#         for line in f_handle:
#             if line[0] == "#":
#                 continue
#             line = line[:-1].split()
#             annotation_type = line[2]
#             if annotation_type == "CDS" and not collecting_exons:
#                 gene_info = " ".join(line[8:]).split(";")[:-1]
#                 info_dic = {e.split()[0]: e.split()[1].replace('"', '')
#                             for e in gene_info}
#                 protein_id = info_dic["protein_id"]
#                 collecting_exons = True
#                 continue


def gene_recall(genome_alignment_dir):
    """Calculate a recall by genes is not possible with bedtools, hence own implementation."""
    gtf_ift = IntervalFrameTree(gtf_file_path=f"{DATA_DIR}/{ANNOTATION_FILE_MOD}", remove_non_coding=True)
    rnacode_ift = IntervalFrameTree(bed_file_path=f"{genome_alignment_dir}/RNAcode.bed")

    hss_score_dic = {}
    for chromosome in os.listdir(genome_alignment_dir):
        if not os.path.isdir(genome_alignment_dir + "/" + chromosome):
            continue
        chromosome_dir_path = f"{genome_alignment_dir}/{chromosome}/"
        hss_score_dic_path = f"{chromosome_dir_path}/hss_score_dic.json"
        with open(hss_score_dic_path, "r", encoding="UTF-8") as f_handle:
            hss_score_dic.update(json.load(f_handle))

    gene_hss_dic = rnacode_ift.overlap_by_gene(gtf_ift)

    p_val_list = []
    not_found_genes = []
    for gene_id, intervals in gene_hss_dic.items():
        if len(intervals) == 0:
            not_found_genes.append(gene_id)
            continue
        # This is the cut off from RNAcode
        best_p = 0.01
        for hss in intervals:
            hss = hss.split()
            hss_id = hss[3]
            # legacy
            # hss_id = f"{hss[0]}_{hss[3]}"

            best_p = min(hss_score_dic[hss_id], best_p)
        p_val_list.append(best_p)

    p_val_list.sort()

    with open(f"{genome_alignment_dir}/gene_recall.tsv", "w", encoding="UTF-8") as f_handle:
        # for i in reversed(np.arange(1, 17, 0.000001)):
        for p_val in np.logspace(-1, -17, num=100):
            # p_val = pow(10, -i)
            f_handle.write(f"{p_val}\t{bisect_right(p_val_list, p_val)}\n")

    with open(f"{genome_alignment_dir}/not_found_genes.lst", "w", encoding="UTF-8") as f_handle:
        f_handle.write("\n".join(not_found_genes) + "\n")


def exon_recall(genome_alignment_dir):
    """Only a redo of bedtools intersect. Same results but my py implementation is much slower."""
    gtf_ift = IntervalFrameTree(gtf_file_path=f"{DATA_DIR}/{ANNOTATION_FILE_MOD}")
    rnacode_ift = IntervalFrameTree(bed_file_path=f"{genome_alignment_dir}/RNAcode.bed")
    overlap_ift = rnacode_ift.overlap(gtf_ift)
    with open(f"{MULTIZ20WAY_DIR}/RNAcode_overlap_PY.bed", "w", encoding="UTF-8") as f_handle:
        f_handle.write(overlap_ift.make_annotation())

    no_overlap_ift = rnacode_ift.no_overlap(gtf_ift)
    with open(f"{MULTIZ20WAY_DIR}/RNAcode_no_overlap_PY.bed", "w", encoding="UTF-8") as f_handle:
        f_handle.write(no_overlap_ift.make_annotation())


def main():
    """Execute the recall calculation."""
    genome_alignment_dir = sys.argv[1]
    if not os.path.isdir(genome_alignment_dir):
        print(f"{genome_alignment_dir} is not a directory")
        sys.exit(1)
    gene_recall(genome_alignment_dir)


if __name__ == "__main__":
    main()

"""

import json
import IntervalFrameTree
from importlib import reload

genome_alignment_dir = "/scr/k61san2/john/rnacode_human_CS/multiz100way/"

genome_alignment_dir = "/scr/k61san2/john/rnacode_human/multiz100way/"

with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR_OLD = parameters["DATA_DIR_OLD"]

DATA_DIR = parameters["DATA_DIR"]

HG38_ANNOTATION_WEB_FTP = parameters["HG38_ANNOTATION_WEB_FTP"]
ANNOTATION_FILE = HG38_ANNOTATION_WEB_FTP.split("/")[-1].replace(".gz", "")
ANNOTATION_FILE_MOD = ANNOTATION_FILE.replace(".gtf", ".cds.gtf")
MULTIZ100WAY_DIR = DATA_DIR + "/multiz100way/"

def info_line_to_dic(info_line):
    info_dic = {}
    for field in info_line.split(";"):
        if field in ["", "\n"]:
            continue
        field = field.split('"')
        info_dic[field[0].replace(" ", "")] = field[1]
    return info_dic

gtf_ift = IntervalFrameTree.IntervalFrameTree(gtf_file_path=f"{DATA_DIR}/{ANNOTATION_FILE_MOD}", remove_non_coding=True)

rnacode_ift = IntervalFrameTree.IntervalFrameTree(bed_file_path=f"{genome_alignment_dir}/RNAcode.bed")

gtf_ift = IntervalFrameTree.IntervalFrameTree(gtf_file_path="./test.gtf", remove_non_coding=True)

gene_dic = {}
for chromosome in rnacode_ift.genome_dic:
    if chromosome not in gtf_ift.genome_dic:
        continue
    for frame in (-3, -2, -1, 1, 2, 3):
        rnacode_ift_it = rnacode_ift.genome_dic[chromosome][frame]
        gtf_ift_it = gtf_ift.genome_dic[chromosome][frame]
        for interval in gtf_ift_it.items():
            line_split = interval.data[1].split("\t")
            info_line = " ".join(line_split[8:])
            info_dic = info_line_to_dic(info_line)
            gene_id = info_dic["gene_id"]
            if gene_id not in gene_dic:
                gene_dic[gene_id] = []
            if gene_id == "ENST00000381317":
                print(frame)
                print(interval[:-1])
            for overlap_interval in rnacode_ift_it.overlap(interval):
                gene_dic[gene_id].append(overlap_interval.data[1])


rnacode_ift = IntervalFrameTree.IntervalFrameTree(bed_file_path=f"{genome_alignment_dir}/RNAcode.bed")

IntervalFrameTree = reload(IntervalFrameTree)

gtf_ift = IntervalFrameTree.IntervalFrameTree(gtf_file_path="./test.gtf", remove_non_coding=True)

for frame in (-3, -2, -1):
    for interval in rnacode_ift.genome_dic["chr7"][frame].items():
        line_split = interval.data[1].split()
        id = line_split[3]
        if id == "HSS_833-hg38.chr7-block_index-1497197":
            hhs_int = interval
            print(frame)
            print(interval[:2])
            break

for frame in (-3, -2, -1):
    for interval in gtf_ift.genome_dic["chr7"][frame].items():
        line_split = interval.data[1].split()
        info_line = " ".join(line_split[8:])
        info_dic = info_line_to_dic(info_line)
        gene_id = info_dic["gene_id"]
        if gene_id == "ENSG00000227191":
            if not info_dic["exon_number"] == "1":
                continue
            print(interval.data[1])
            gen_int = interval
            # print(info_dic["exon_number"])
            print(frame)
            print(interval[:2])

gen_int

for overlap_interval in self_it.overlap(interval):
    print(overlap_interval)
pos = (18773126, 18773404)
print((pos[1] - pos[0]) % 3)

pos = (6868036, 6868462)
print((pos[1] - pos[0]) % 3)


6868203 % 3

6868037 % 3

-1
Interval(6868203, 6868263, (4620, 'chrY 6868203 6868263 HSS_1365-chrY_298_22 7 - 6868203 6868263 255,128,

"""
