#!/usr/bin/python3
"""Find not overlapping HSS to inspect."""

import json
import os
from glob import glob
import subprocess


with open("../parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
DATA_DIR_OLD = parameters["DATA_DIR_OLD"]
genome_alignment_dir = DATA_DIR + "/multiz100way/"
false_anno_path = genome_alignment_dir + "/RNAcode_no_overlap.bed"

# only the last part needs to be changed if a new UCSC session was used
base_url = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={}%3A{}%2D{}&hgsid=289764578_6cpu3ma1XGs6ShjdAiG0hltYiGwb"

hss_score_dic = {}
for chromosome in os.listdir(genome_alignment_dir):
    if not os.path.isdir(genome_alignment_dir + "/" + chromosome):
        continue
    chromosome_dir_path = f"{genome_alignment_dir}/{chromosome}/"
    hss_score_dic_path = f"{chromosome_dir_path}/hss_score_dic.json"
    with open(hss_score_dic_path, "r", encoding="UTF-8") as f_handle:
        hss_score_dic.update(json.load(f_handle))

i = 0
with open(false_anno_path, "r", encoding="UTF-8") as f_handle:
    for line in f_handle:
        if line[0] == "#":
            continue
        if line[0] == "\n":
            continue

        line_split = line.split()
        chromosome = line_split[0]
        start = int(line_split[1])
        end = int(line_split[2])
        hss_id = line_split[3]

        if float(hss_score_dic[hss_id]) > 0.000000001:
            continue

        if chromosome == "chrY":
            continue
        print(i)
        i += 1
        if i < 13:
            continue
        if i > 30:
            break

        # if "split" in hss_id:
        #     continue
        # if "," in hss_id:
        #     continue

        block_id = "-".join(hss_id.split("-")[1:])

        found = False
        maf_block = "a scrore=0"
        for maf_file in glob(f"{genome_alignment_dir}/{chromosome}/big_blocks/*maf"):
            with open(maf_file, "r", encoding="UTF-8") as f_handle:
                for line in f_handle:
                    if block_id in line:
                        maf_block += line
                        found = True
                    if found:
                        maf_block += line
                        if line == "\n":
                            break
            if found:
                break

        with open("./test.maf", "w", encoding="UTF-8") as f_handle:
            f_handle.write(maf_block)

        call_str = "RNAcode -e test.maf"
        try:
            completed_process = subprocess.run(
                call_str.split(), capture_output=True, text=True, check=True
            )
        except subprocess.CalledProcessError as exc:
            print(exc.stdout)
            print(exc.stderr)
            raise

        print(hss_id)
        print(hss_score_dic[hss_id])
        print(int(end)-int(start))
        print(base_url.format(chromosome, start, end))
        print()
        call_str = "evince eps/hss-0.eps"
        try:
            completed_process = subprocess.run(
                call_str.split(), capture_output=True, text=True, check=True
            )
        except subprocess.CalledProcessError as exc:
            print(exc.stdout)
            print(exc.stderr)
            raise
