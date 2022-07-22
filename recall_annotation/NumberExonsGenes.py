#!/usr/bin/python3
"""Caculate the number of exons and gene in the annotation.

Needed for plottin in Rscript (./fdr_recal.R) the number is simply copy pasted
into the script.
"""

import json


with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
HG38_ANNOTATION_WEB_FTP = parameters["HG38_ANNOTATION_WEB_FTP"]

ANNOTATION_FILE = HG38_ANNOTATION_WEB_FTP.split("/")[-1].replace(".gz", "")
ANNOTATION_FILE_MOD = ANNOTATION_FILE.replace(".gff", ".mod.gff")

exon_set = set()
gene_set = set()
with open(f"{DATA_DIR}/{ANNOTATION_FILE_MOD}", "r", encoding="UTF-8") as f_handle:
    for line in f_handle:
        if line[0] == "#":
            continue
        line_split = line.split()
        annotation_type = line_split[2]
        gene_info = " ".join(line_split[8:]).split(";")[:-1]
        info_dic = {e.split()[0]: e.split()[1].replace('"', '')
                    for e in gene_info}
        if annotation_type == "exon":
            exon_set.add(info_dic["exon_id"])
            gene_set.add(info_dic["gene_id"])
print(f"Number of exons {len(exon_set)}")
print(f"Number of genes {len(gene_set)}")
