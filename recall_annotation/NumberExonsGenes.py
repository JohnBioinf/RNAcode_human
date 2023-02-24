#!/usr/bin/python3
"""Caculate the number of exons and gene in the annotation.

Needed for plottin in Rscript (./fdr_recal.R) the number is simply copy pasted
into the script.
"""

import json
from IntervalFrameTree import ENSEMBLE_CODING_BIOTYPES
from IntervalFrameTree import ENSEMBLE_PSEUDO_CODING_BIOTYPES
from IntervalFrameTree import info_line_to_dic

with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
HG38_ANNOTATION_WEB_FTP = parameters["HG38_ANNOTATION_WEB_FTP"]

ANNOTATION_FILE = HG38_ANNOTATION_WEB_FTP.split("/")[-1].replace(".gz", "")
ANNOTATION_FILE_CDS = ANNOTATION_FILE.replace(".gtf", ".cds.gtf")
ANNOTATION_FILE_PSEUDO = ANNOTATION_FILE.replace(".gtf", ".pseudo.gtf")

gene_set = set()
# biotype_set = set()
out_pseudo_handle = open(f"{DATA_DIR}/{ANNOTATION_FILE_PSEUDO}", "w", encoding="UTF-8")
out_cds_handle = open(f"{DATA_DIR}/{ANNOTATION_FILE_CDS}", "w", encoding="UTF-8")
with open(f"{DATA_DIR}/{ANNOTATION_FILE}", "r", encoding="UTF-8") as f_handle:
    for line in f_handle:
        if line[0] == "#":
            out_pseudo_handle.write(line)
            out_cds_handle.write(line)
            continue
        line_split = line.split()
        annotation_type = line_split[2]

        info_line = " ".join(line_split[8:])
        gene_info_dic = info_line_to_dic(info_line)
        # biotype_set.add(gene_info_dic["gene_biotype"])
        if "gene_biotype" in gene_info_dic:
            if gene_info_dic["gene_biotype"] in ENSEMBLE_PSEUDO_CODING_BIOTYPES:
                out_pseudo_handle.write("chr" + line)
                continue
        else:
            print("GTF line has not biotype. Can not discriminate if coding or not")
            print(line)
            exit(1)

        if annotation_type == "CDS" and gene_info_dic["transcript_biotype"] in ENSEMBLE_CODING_BIOTYPES:
            gene_set.add(gene_info_dic["gene_id"])
            out_cds_handle.write("chr" + line)

# print(list(biotype_set))

out_cds_handle.close()
out_pseudo_handle.close()

print(f"Number of genes {len(gene_set)}")
