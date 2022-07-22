#!/usr/bin/python3
"""Calculates the recall rate and FDR for RNAcode annotation.

The calculaten works under the assumption that the current annoation is the ground truth.
"""
import json
from IntervalFrameTree import IntervalFrameTree
import os
from bisect import bisect_right
import sys


with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
HG38_ANNOTATION_WEB_FTP = parameters["HG38_ANNOTATION_WEB_FTP"]

ANNOTATION_FILE = HG38_ANNOTATION_WEB_FTP.split("/")[-1].replace(".gz", "")
ANNOTATION_FILE_MOD = ANNOTATION_FILE.replace(".gff", ".mod.gff")
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
    gtf_ift = IntervalFrameTree(gtf_file_path=f"{DATA_DIR}/{ANNOTATION_FILE_MOD}")
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
    for _gene_id, intervals in gene_hss_dic.items():
        if len(intervals) == 0:
            continue
        best_p = 0.01
        for hss in intervals:
            hss = hss.split()
            best_p = min(hss_score_dic[f"{hss[0]}_{hss[3]}"], best_p)
        p_val_list.append(best_p)

    uniq_p_val = list(set(p_val_list))
    uniq_p_val.sort()
    p_val_list.sort()

    with open(f"{genome_alignment_dir}/gene_recall.tsv", "w", encoding="UTF-8") as f_handle:
        for p_val in uniq_p_val:
            f_handle.write(f"{p_val}\t{bisect_right(p_val_list, p_val)}\n")


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
