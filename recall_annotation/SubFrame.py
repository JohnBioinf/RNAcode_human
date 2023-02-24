#!/usr/bin/python3
"""Generates a Subset of an annotation file based only containing entrys in a specific frame.

When the input is a GTF file also only exons will be selected.
"""

import sys
from IntervalFrameTree import IntervalFrameTree


def main():
    """Get answer of Life, the Universe and Everything."""
    frame = int(sys.argv[1])
    annotation_file = sys.argv[2]
    if ".gtf" in annotation_file:
        IntervalFrameTree(gtf_file_path=annotation_file, sub_frame=frame)
    elif ".bed" in annotation_file:
        IntervalFrameTree(bed_file_path=annotation_file, sub_frame=frame)
    else:
        print("Unknown file type")
        sys.exit(1)


if __name__ == "__main__":
    main()
