#!/usr/bin/python3
"""Test intervalFrameTree."""

from IntervalFrameTree import IntervalFrameTree


def main():
    """Test intervalFrameTree."""
    gtf_file_path = "./test.gtf"
    gtf_ift = IntervalFrameTree(gtf_file_path=gtf_file_path, remove_non_coding=True)
    for _chrom, frame_dic in gtf_ift.genome_dic.items():
        for frame, IT in frame_dic.items():
            for interval in IT:
                print(frame)
                print(interval)
                print()


if __name__ == "__main__":
    main()
