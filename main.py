#!/usr/bin/python3
"""Pipeline that computes RNAcode on a genome wide alignment in parallel.

The pipeline will first concat maf blocks, for which a concatination is trivial
and the distances are not to great (12nucs/ 4 AS). These larger blocks will than
be splited with a sliding window and be packed in new maf file (big blocks).
After this an anlysis with RNAcode is performed. Failed big blocks are repeated
and if they still failled, split up into single maf blocks per file. And
anlaysed individualy. All maf blocks which failed twice indvidually will be
reported in the file './bad_maf_blocks.txt'.
The results are than used to build a bed file.

The concatination and RNAcode analysis are performed in parallel.
"""

import json
import os
import sys
import subprocess
from glob import glob
from datetime import datetime
import re
from math import log

from MafBlock import MafBlock


with open("./parameters_local.json", "r", encoding="UTF-8") as file_handle:
    parameters = json.load(file_handle)

DATA_DIR = parameters["DATA_DIR"]
DATA_DIR_OLD = parameters["DATA_DIR_OLD"]
MULTIZ100WAY_WEB_FTP = parameters["MULTIZ100WAY_WEB_FTP"]
MULTIZ20WAY_WEB_FTP = parameters["MULTIZ20WAY_WEB_FTP"]

# Parameters for preprocessing can be found stream_chromosome.py
P_THRESHOLD = 0.01
NUM_CPUS = 20

MULTIZ100WAY_DIR = DATA_DIR + "/multiz100way/"
MULTIZ100WAY_DIR_OLD = DATA_DIR_OLD + "/multiz100way/"
MULTIZ20WAY_DIR = DATA_DIR + "/multiz20way/"


def eprint(*a, **k):
    """Print to stderr."""
    print(*a, file=sys.stderr, **k)


def system_call(call_str, env=None):
    """Perform system call based on input str.

    :param str call_str: The command as a string.
    :return: The return code and the output
    :rtype: tuple
    """
    with subprocess.Popen(
        call_str.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env
    ) as process:
        process.wait()
        return process.returncode, process.communicate()[0].decode("UTF-8")


def download_maf_file(chromosome_dir, maf_file, ftp_url, checksum_maf_file):
    """Download maf file."""
    if os.path.isfile(f"{chromosome_dir}/{maf_file}"):
        call = f"md5sum {chromosome_dir}/{maf_file}"
        md5sum_maf = system_call(call)[1].split()[0]
        if md5sum_maf != checksum_maf_file:
            print(f"{maf_file} corupt! Redownload")
            os.remove(f"{chromosome_dir}/{maf_file}")
            out = system_call(f"wget -nv -P {chromosome_dir} {ftp_url}/{maf_file}")
            if out[0] != 0:
                eprint("ERROR!")
                eprint(f"Downloading {maf_file} failed")
                eprint(out[1])
                sys.exit(1)
            else:
                print(f"{maf_file} downloaded.")
            call = f"md5sum {chromosome_dir}/{maf_file}"
            md5sum_maf = system_call(call)[1].split()[0]
            if md5sum_maf != checksum_maf_file:
                print(f"{maf_file} corupt!")
                sys.exit(1)
        else:
            print(f"{maf_file} already downloaded")
    else:
        print(f"{maf_file} not downloaded")
        out = system_call(f"wget -nv -P {chromosome_dir} {ftp_url}/{maf_file}")
        if out[0] != 0:
            eprint("ERROR!")
            eprint(f"Downloading {maf_file} failed")
            eprint(out[1])
            sys.exit(1)
        else:
            print(out[1])
            print(f"{maf_file} downloaded.")

        call = f"md5sum {chromosome_dir}/{maf_file}"
        md5sum_maf = system_call(call)[1].split()[0]
        if md5sum_maf != checksum_maf_file:
            print(f"{maf_file} corupt!")
            sys.exit(1)
    sys.stdout.flush()


def init_work_dir(dir_path, ftp_url):
    """Set up the working directory."""
    print("Init working directory.")
    check_sum_file_path = f"{dir_path}/md5sum.txt"

    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)

    if not os.path.isfile(check_sum_file_path):
        call = f"wget -nv -P {dir_path} {ftp_url}/md5sum.txt"
        out = system_call(f"wget -nv -P {dir_path} {ftp_url}/md5sum.txt")
        if out[0] != 0:
            eprint("ERROR!")
            eprint("Downloading chechsum file failed")
            eprint(out[1])
            eprint(call)
            sys.exit(1)

    checksum_dic = {}
    with open(check_sum_file_path, "r", encoding="UTF-8") as f_handle:
        checksum_dic = {line.split()[1]: line.split()[0] for line in f_handle}

    out = system_call(f"curl -ls {ftp_url}")
    if out[0] != 0:
        eprint("ERROR!")
        eprint("Retrieving file list failed!")
        eprint(out[1])
        sys.exit(1)
    maf_file_list = out[1].split("\n")

    for maf_file in maf_file_list:
        if ".maf.gz" not in maf_file:
            continue
        chromosome = maf_file.split(".")[0]
        chromosome_dir = f"{dir_path}/{chromosome}"
        if not os.path.isdir(chromosome_dir):
            os.mkdir(chromosome_dir)
        download_maf_file(chromosome_dir, maf_file, ftp_url, checksum_dic[maf_file])
    print("Finished building working directory")


def compute_genome_alignment_big_blocks(genome_alignment_dir):
    """Compute maf file."""
    print("Concating and spliting genome.")
    start_time_maf_stream = datetime.now()
    call_str = f"./stream_chromosome_parallel.sh {genome_alignment_dir} {NUM_CPUS}"
    completed_process = subprocess.run(call_str.split(), capture_output=True, text=True)
    if completed_process.returncode != 0:
        eprint(completed_process.stderr)
        raise SystemError("Maf stream parallel failed!")
    print(f"Finished preprocessing after {datetime.now() - start_time_maf_stream}.")

    start_time_rnacode = datetime.now()
    call_str = f"./RNAcode_parallel.sh {genome_alignment_dir} big_blocks {NUM_CPUS}"
    completed_process = subprocess.run(call_str.split(), capture_output=True, text=True)
    if completed_process.returncode != 0:
        eprint(completed_process.stderr)
        raise SystemError("RNAcode parallel failed!")

    print(
        f"Finished RNAcode analysis after {datetime.now() - start_time_rnacode}."
    )


def check_errors_rnacode(genome_alignment_dir):
    """Loop over error messages from rnacode and checks for new ones."""
    for chromosome in os.listdir(genome_alignment_dir):
        if not os.path.isdir(genome_alignment_dir + "/" + chromosome):
            continue
        chromosome_dir_path = f"{genome_alignment_dir}/{chromosome}/"
        for err_file in glob(f"{chromosome_dir_path}/big_blocks/*err"):
            line_before = "First line"
            with open(err_file, "r", encoding="UTF-8") as f_handle:
                for line in f_handle:
                    if line == "\n":
                        continue
                    if "Target: " in line:
                        line_before = line
                        continue
                    if "ERROR: Empty alignment file" in line:
                        line_before = line
                        continue
                    print(chromosome)
                    print(err_file)
                    print(line_before)
                    print(line)


def find_failed_twice(chromosome_dir_path, block_type="big_block"):
    """Check parallel log for double fails."""
    if block_type == "big_block":
        reg_pat = re.compile(r"big_block_[0-9]+\.maf")
        log_file = f"{chromosome_dir_path}/big_blocks_parallel.log"
    elif block_type == "single_block":
        reg_pat = re.compile(r"big_block_[0-9]+-s_[0-9]+\.maf")
        log_file = f"{chromosome_dir_path}/single_blocks_parallel.log"
    else:
        raise ValueError(
            f'Block type must be either "big_block" or "single_block". Not {block_type}'
        )

    failed_blocks = []
    failed_twice = []
    if not os.path.isfile(log_file):
        return failed_twice

    with open(log_file, "r", encoding="UTF-8") as f_handle:
        log_file_content = f_handle.read().split("\n")
        for line in log_file_content[1:-1]:
            exitval = int(line.split()[6])
            if exitval != 0:
                try:
                    block = reg_pat.findall(line)[0].split(".")[0]
                except IndexError:
                    print(f"Error can not find {block_type.replace('_', ' ')} in line")
                    print(line)
                    continue
                    # sys.exit(1)
                if block in failed_blocks:
                    failed_twice.append(block)
                else:
                    failed_blocks.append(block)

    return failed_twice


def build_single_blocks(chromosome_dir_path, failed_twice):
    """Split big blocks that failed into single maf-blocks."""
    block_index = 0
    maf = MafBlock()
    for big_block in failed_twice:
        with open(
            f"{chromosome_dir_path}/big_blocks/{big_block}.maf",
            "r",
            encoding="UTF-8",
        ) as f_handle:
            for line in f_handle:
                if line == "\n":
                    block_index += 1
                    file_path = f"{chromosome_dir_path}/big_blocks/{big_block}-s_{block_index}.maf"

                    with open(file_path, "w", encoding="UTF-8") as f_handle:
                        f_handle.write(str(maf))
                    maf = MafBlock()
                else:
                    maf.add(line)


def check_failed_and_retry(genome_alignment_dir):
    """Check logs of all big maf-blocks and retry as single block."""
    print("Check for big blocks which failed twice")
    start_time_rnacode = datetime.now()
    for chromosome in os.listdir(genome_alignment_dir):
        if not os.path.isdir(genome_alignment_dir + "/" + chromosome):
            continue

        print(chromosome, flush=True)
        start_time_chromosome = datetime.now()
        chromosome_dir_path = f"{genome_alignment_dir}/{chromosome}/"

        failed_twice = find_failed_twice(chromosome_dir_path)
        if len(failed_twice) == 0:
            print("No big block failed twice")
            continue

        print("The following blocks failed twice:")
        print("\n".join(failed_twice))

        build_single_blocks(chromosome_dir_path, failed_twice)
        print("Compute single blocks.")
        call = f"./RNAcode_parallel.sh {chromosome_dir_path} single_blocks {NUM_CPUS}"
        out = system_call(call)
        if out[0] != 0:
            eprint("Error parallel")
            eprint(out[1])
        print(
            f"Finished chromosome single blocks after {datetime.now() - start_time_chromosome}."
        )
    print(
        f"Finished computing genome allignment after {datetime.now() - start_time_rnacode}."
    )


def build_segments(rnacode_res, chromosome, only_best=True):
    """Build segments from rnacode results as intermediate structure to bed line."""
    segments = []
    for line in rnacode_res:
        start = int(line[7])
        end = int(line[8])
        hss_id = line[0]
        maf_block = line[6]
        strand = line[1]
        p_val = float(line[-1])

        new_segment = {
            "start": start,
            "end": end,
            "strand": strand,
            "id": f"HSS_{hss_id}-{maf_block}",
            "chromosome": chromosome,
            "p_val": p_val,
        }
        for segment in segments:
            # Check overlap
            if (
                new_segment["start"] <= segment[0]["end"]
                and segment[0]["start"] <= new_segment["end"]
            ):
                # Not same frame
                if not (
                    segment[0]["strand"] == new_segment["strand"]
                    and segment[0]["start"] % 3 == new_segment["start"] % 3
                ):
                    # if only best should be kept check if new segment is better
                    # than old one
                    if only_best and segment[0]["p_val"] > new_segment["p_val"]:
                        segment[0] = new_segment
                        break
                # Same frame
                else:
                    segment[0]["start"] = min((segment[0]["start"], new_segment["start"]))
                    segment[0]["end"] = max((segment[0]["end"], new_segment["end"]))
                    segment[0]["id"] += "," + new_segment["id"]
                    segment[0]["p_val"] = min(new_segment["p_val"], segment[0]["p_val"])
                    break
        # if not overlap found simply add
        else:
            segments.append([new_segment])
    return [seg[0] for seg in segments]



def build_bed(genome_alignment_dir):
    """Build bed file from RNAcode results."""
    print("Build bed file")
    genome_bed_file_path = f"{genome_alignment_dir}/RNAcode.bed"
    if os.path.isfile(genome_bed_file_path):
        os.remove(genome_bed_file_path)
    bad_maf_file_path = "./bad_maf_blocks.txt"
    if os.path.isfile(bad_maf_file_path):
        os.remove(bad_maf_file_path)

    for chromosome in os.listdir(genome_alignment_dir):
        if not os.path.isdir(genome_alignment_dir + "/" + chromosome):
            continue
        if chromosome != "chr1":
            continue
        print(f"Processing {chromosome}")
        chromosome_dir_path = f"{genome_alignment_dir}/{chromosome}/"
        bed_file_path = f"{chromosome_dir_path}/RNAcode.bed"
        hss_score_dic_path = f"{chromosome_dir_path}/hss_score_dic.json"

        # failed_twice = find_failed_twice(chromosome_dir_path, block_type="single_block")
        # # Note single maf blocks which failed twice
        # with open(bad_maf_file_path, "a", encoding="UTF-8") as f_handle:
        #     f_handle.write(
        #         "\n".join(
        #             [
        #                 f"{genome_alignment_dir}big_blocks/{mb}.maf"
        #                 for mb in failed_twice
        #             ]
        #         )
        #         + "\n"
        #     )
        # failed_twice += find_failed_twice(chromosome_dir_path, block_type="big_block")
        # # if len(failed_twice) != 0:
        # #     print(f"The following blocks failed twice {', '.join(failed_twice)}")

        print("Collect RNAcode results")
        rnacode_res = []
        for rnacode_res_file_path in glob(f"{chromosome_dir_path}/big_blocks/*res.tsv"):
            # block = rnacode_res_file_path.split("/")[-1].replace(".maf.res.tsv", "")
            # if block in failed_twice:
            #     # print(f"{block} failed twice.")
            #     continue
            with open(rnacode_res_file_path, "r", encoding="UTF-8") as f_handle:
                for line in f_handle.read().split("\n")[:-1]:
                    line = line.split("\t")
                    if float(line[-1]) > P_THRESHOLD:
                        continue
                    rnacode_res.append(line)

        print("Build segments")
        segments = build_segments(rnacode_res, chromosome)

        hss_score_dic = {
            seg["id"]: seg["p_val"]
            for seg in segments
        }

        with open(hss_score_dic_path, "w", encoding="UTF-8") as f_handle:
            json.dump(hss_score_dic, f_handle)

        print("Make bed lines")
        bed_lines = hss_to_bed_line(segments)

        with open(bed_file_path, "w", encoding="UTF-8") as f_handle:
            f_handle.write(
                "\n".join([" ".join(map(str, line)) for line in bed_lines]) + "\n"
            )

        with open(genome_bed_file_path, "a", encoding="UTF-8") as f_handle:
            f_handle.write(
                "\n".join([" ".join(map(str, line)) for line in bed_lines]) + "\n"
            )


def hss_to_bed_line(segments):
    """Convert a list of high scoring segments from RNAcode to a bed line."""
    block_count = 1
    block_starts = 0
    exp_ids = 1
    exp_count = 1

    frame_color_dic = {}
    # frame 1 red
    frame_color_dic[1] = ["255,128,128", "255,26,26"]
    # frame 2 blue
    frame_color_dic[2] = ["128,128,255", "26,26,255"]
    # frame 3 green
    frame_color_dic[3] = ["153,230,153", "45,185,45"]
    # frame 4 purple
    frame_color_dic[4] = ["223,128,255", "172,0,230"]
    # frame 5 orange
    frame_color_dic[5] = ["255,191,128", "230,115,0"]
    # frame 6 turquoise
    frame_color_dic[6] = ["128,255,229", "0,230,184"]

    lines = []
    for segment in segments:
        p_val = segment["p_val"] if segment["p_val"] > 0 else 0.000001
        score = min(int(log(p_val, 2) * -1), 1000)
        name = segment["id"]
        thick_end = chrom_end = segment["end"]
        thick_start = chrom_start = segment["start"]
        if segment["strand"] == "-1":
            frame = chrom_start % 3 + 4
        else:
            frame = chrom_start % 3 + 1
        quality = 0 if segment["p_val"] > 0.001 else 1

        item_rgb = frame_color_dic[frame][quality]
        block_sizes = chrom_end - chrom_start
        lines.append(
            [
                segment["chromosome"],
                chrom_start,
                chrom_end,
                name,
                score,
                segment["strand"],
                thick_start,
                thick_end,
                item_rgb,
                block_count,
                block_sizes,
                block_starts,
                exp_count,
                exp_ids,
            ]
        )
    return lines


def full_pipeline(work_dir, web_ftp):
    """Full Piepline."""
    # init_work_dir(work_dir, web_ftp)
    # compute_genome_alignment_big_blocks(work_dir)
    # check_failed_and_retry(work_dir)
    build_bed(work_dir)


def main():
    """Execute pipeline according to global parameters."""
    full_pipeline(MULTIZ100WAY_DIR, MULTIZ100WAY_WEB_FTP)
    full_pipeline(MULTIZ100WAY_DIR_OLD, MULTIZ100WAY_WEB_FTP)


if __name__ == "__main__":
    main()
