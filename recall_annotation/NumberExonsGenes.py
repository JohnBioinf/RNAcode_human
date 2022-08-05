#!/usr/bin/python3
"""Caculate the number of exons and gene in the annotation.

Needed for plottin in Rscript (./fdr_recal.R) the number is simply copy pasted
into the script.
"""

import json
from IntervalFrameTree import ENSEMBLE_CODING_BIOTYPES
from IntervalFrameTree import info_line_to_dic


with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
HG38_ANNOTATION_WEB_FTP = parameters["HG38_ANNOTATION_WEB_FTP"]

ANNOTATION_FILE = HG38_ANNOTATION_WEB_FTP.split("/")[-1].replace(".gz", "")
ANNOTATION_FILE_MOD = ANNOTATION_FILE.replace(".gtf", ".mod.gtf")

exon_count = 0
gene_set = set()
out_f_handle = open(f"{DATA_DIR}/{ANNOTATION_FILE_MOD}", "w", encoding="UTF-8")
with open(f"{DATA_DIR}/{ANNOTATION_FILE}", "r", encoding="UTF-8") as f_handle:
    for line in f_handle:
        if line[0] == "#":
            out_f_handle.write(line)
            continue
        line_split = line.split()
        annotation_type = line_split[2]

        if annotation_type == "CDS":
            info_line = " ".join(line_split[8:])
            gene_info_dic = info_line_to_dic(info_line)
            if "transcript_biotype" in gene_info_dic:
                if gene_info_dic["transcript_biotype"] not in ENSEMBLE_CODING_BIOTYPES:
                    continue
            else:
                print("GTF line has not biotype. Can not discriminate if coding or not")
                print(line)
                exit(1)

            out_f_handle.write("chr" + line)

            exon_count += 1
            gene_set.add(gene_info_dic["gene_id"])

out_f_handle.close()

print(f"Number of exons {exon_count}")
print(f"Number of genes {len(gene_set)}")
