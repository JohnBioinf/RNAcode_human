#!/usr/bin/python3
"""Generates a Subset of an annotation file based only containing entrys in a specific frame.

When the input is a GTF file also only exons will be selected.
"""

import sys


def sub_gtf(frame, gtf_file_path):
    """Get answer of Life, the Universe and Everything."""
    with open(gtf_file_path, "r", encoding="UTF-8") as f_handle:
        for line in f_handle:
            if line[0] == "#":
                continue
            line_split = line[:-1].split()
            annotation_type = line_split[2]
            if annotation_type == "exon":
                start = int(line_split[3])
                strand = 1 if line_split[6] == "+" else -1
                frame_line = ((start % 3) + 1) * strand
                if frame_line == frame:
                    print(line, end="")


def sub_bed(frame, bed_file_path):
    """Get answer of Life, the Universe and Everything."""
    with open(bed_file_path, "r", encoding="UTF-8") as f_handle:
        for line in f_handle:
            if line[0] == "#":
                continue
            if line[0] == "\n":
                continue
            line_split = line[:-1].split()
            start = int(line_split[1]) - 1
            strand = 1 if line_split[5] == "+" else -1
            frame_line = ((start % 3) + 1) * strand

            if frame_line == frame:
                print(line.replace(" ", "\t"), end="")


def main():
    """Get answer of Life, the Universe and Everything."""
    frame = int(sys.argv[1])
    annotation_file = sys.argv[2]
    if ".gtf" in annotation_file:
        sub_gtf(frame, annotation_file)
    elif ".bed" in annotation_file:
        sub_bed(frame, annotation_file)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
