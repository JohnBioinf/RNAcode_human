#!/usr/bin/python3
"""Python module for the class IntervalFrameTree.

The IntervallFrameTree provides a class that can be used to find overlapping
genomic regions based on an interval tree. In contrast to a normal interval tree
it is frame and chromosome sensitiv, hence two intervals are only overlapping if
they are on the same chromosome and frame.
"""
from intervaltree import IntervalTree  # , Interval
import sys

ENSEMBLE_CODING_BIOTYPES = ["protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TEC", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
ENSEMBLE_PSEUDO_CODING_BIOTYPES = ['transcribed_unprocessed_pseudogene', 'transcribed_processed_pseudogene', 'unprocessed_pseudogene', 'pseudogene', 'polymorphic_pseudogene', 'TR_V_pseudogene', 'TEC', 'translated_processed_pseudogene', 'IG_V_pseudogene', 'IG_pseudogene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'processed_pseudogene', 'IG_C_pseudogene', 'transcribed_unitary_pseudogene', 'unitary_pseudogene', 'translated_unprocessed_pseudogene']


def info_line_to_dic(info_line):
    """Build a dictinoray from the info line of a gtf file."""
    info_dic = {}
    for field in info_line.split(";"):
        if field in ["", "\n"]:
            continue
        field = field.split('"')
        info_dic[field[0].replace(" ", "")] = field[1]
    return info_dic


class IntervalFrameTree:
    """An interval tree for a genome, called interval frame tree (IFT).

    The IntervallFrameTree provides a class that can be used to find overlapping
    genomic regions based on an interval tree. In contrast to a normal interval tree
    it is frame and chromosome sensitiv, hence two intervals are only overlapping if
    they are on the same chromosome and frame.
    """

    def __init__(self, gtf_file_path=None, bed_file_path=None, annotation_type=None, remove_non_coding=False, sub_frame=None):
        """Init an empty IFT.

        Usage is to give the path to the annotation file. Or only give the
        annotation type to make an empty IFT.
        Parameter remove_non_coding is "True" remove all non-coding entrys.
        Parameter sub_frame is set do not build IFT only print lines in the
        given frame.
        """
        self.genome_dic = {}
        self.comment_lines = []
        self.annotation_type = None
        if sub_frame not in [None, 1, 2, 3, -1, -2, -3]:
            raise ValueError("The sub frame parameter must be between -3 and 3.")
        if gtf_file_path:
            self.annotation_type = "gtf"
            if self.annotation_type != annotation_type and annotation_type:
                raise ValueError(f"The annotation type {annotation_type} doe not match {gtf_file_path}.")
            self.__add_gtf_file(gtf_file_path, sub_frame=sub_frame)
        elif bed_file_path:
            self.annotation_type = "bed"
            if self.annotation_type != annotation_type and annotation_type:
                raise ValueError(f"The annotation type {annotation_type} doe not match {bed_file_path}.")
            with open(bed_file_path, "r", encoding="UTF-8") as f_handle:
                for i, line in enumerate(f_handle):
                    if line[0] == "#":
                        self.comment_lines.append((i, line))
                        continue
                    if line[0] == "\n":
                        continue
                    self.__add_bed_line(i, line, sub_frame=sub_frame)
        else:
            self.annotation_type = annotation_type

    def __add_gtf_file(self, gtf_file_path, remove_non_coding=True, sub_frame=None):
        current_protein_id = None
        offset = 0
        with open(gtf_file_path, "r", encoding="UTF-8") as f_handle:
            for i, line in enumerate(f_handle):
                if line[0] == "#":
                    self.comment_lines.append((i, line))
                    continue
                line_split = line[:-1].split()
                annotation_type = line_split[2]
                if annotation_type == "CDS":
                    info_line = " ".join(line_split[8:])
                    gene_info_dic = info_line_to_dic(info_line)
                    if "transcript_biotype" in gene_info_dic:
                        if gene_info_dic["transcript_biotype"] not in ENSEMBLE_CODING_BIOTYPES:
                            break
                    else:
                        print("GTF line has not biotype. Can not discriminate if coding or not")
                        print(line)
                        sys.exit(1)
                else:
                    continue
                chromosome = line_split[0]
                # The first base is based 1 in gtf the interval tree first base is 0
                start = int(line_split[3]) - 1
                end = int(line_split[4])
                strand = 1 if line_split[6] == "+" else -1
                # The offset is used to find the start of the first coding nucleotide.
                if current_protein_id != gene_info_dic["protein_id"]:
                    # For proteins on positive strand no offset for first exon
                    # to define the current frame.
                    if strand == 1:
                        offset = 0
                    else:
                        offset = (end - start) % 3
                    current_protein_id = gene_info_dic["protein_id"]
                    last_exon_num = int(gene_info_dic["exon_number"])
                else:
                    if last_exon_num > int(gene_info_dic["exon_number"]):
                        print("The gtf file is not ordered according by exons!")
                        sys.exit(1)
                    last_exon_num = int(gene_info_dic["exon_number"])
                    if strand == 1:
                        offset = (end - start - offset) % 3
                    else:
                        offset = (end + offset - start) % 3

                if strand == 1:
                    frame = (((start - offset) % 3) + 1) * strand
                else:
                    frame = (((start + offset) % 3) + 1) * strand

                # For proteins on positive strand offset for second exon.
                if last_exon_num == 1 and strand == 1:
                    offset = (end - start) % 3

                if sub_frame:
                    if frame == sub_frame:
                        print(line, end="")
                else:
                    if chromosome not in self.genome_dic:
                        self.__add_chromosome(chromosome)
                    # TODO
                    # # Some exon seem to have start and end exactly the same. No clue why.
                    # if start < end:
                    #     self.genome_dic[chromosome][frame][start:end] = (i, line)
                    self.genome_dic[chromosome][frame][start:end] = (i, line)

    def __add_bed_line(self, i,  line, sub_frame):
        """Add bed line to IFT."""
        line_split = line[:-1].split()
        chromosome = line_split[0]
        start = int(line_split[1])
        end = int(line_split[2])
        strand = 1 if line_split[5] == "+" else -1
        frame = ((start % 3) + 1) * strand
        if sub_frame:
            if frame == sub_frame:
                print(line.replace(" ", "\t"), end="")
        else:
            if chromosome not in self.genome_dic:
                self.__add_chromosome(chromosome)
            # TODO
            # if start < end:
            #     self.genome_dic[chromosome][frame][start:end] = (i, line)
            self.genome_dic[chromosome][frame][start:end] = (i, line)

    def __add_chromosome(self, chromosome):
        self.genome_dic[chromosome] = {k: IntervalTree() for k in (-3, -2, -1, 1, 2, 3)}

    def add_interval_tree(self, chromosome, frame, interval_tree):
        """Add interval to IFT."""
        if chromosome not in self.genome_dic:
            self.__add_chromosome(chromosome)
        self.genome_dic[chromosome][frame] = interval_tree

    def __str__(self):
        """Represent class as string."""
        if self.genome_dic == {}:
            return "Empty interval frame tree."
        return f"Interval frame tree with {len(self.genome_dic)} chromosomes."

    def overlap(self, ift):
        """Calculate the overlap between self and another interval tree."""
        if not isinstance(ift, IntervalFrameTree):
            raise TypeError("Input needs to be an interval frame tree.")

        if set(self.genome_dic.keys()).intersection(set(ift.genome_dic.keys())) == 0:
            raise KeyError("Both trees do not share any chromosomes.")

        overlap_ift = IntervalFrameTree(annotation_type=self.annotation_type)
        for chromosome in self.genome_dic:
            if chromosome not in ift.genome_dic:
                continue
            for frame in (-3, -2, -1, 1, 2, 3):
                self_it = self.genome_dic[chromosome][frame]
                other_it = ift.genome_dic[chromosome][frame]
                overlap_it = IntervalTree()
                for interval in other_it.items():
                    for overlap_interval in self_it.overlap(interval):
                        overlap_it.add(overlap_interval)
                overlap_ift.add_interval_tree(chromosome, frame, overlap_it)
        return overlap_ift

    def no_overlap(self, other):
        """Calculate the overlap between self and another interval tree."""
        if not isinstance(other, IntervalFrameTree):
            raise TypeError("Input needs to be an interval frame tree.")

        if set(self.genome_dic.keys()).intersection(set(other.genome_dic.keys())) == 0:
            raise KeyError("Both trees do not share any chromosomes.")

        no_overlap_ift = self.overlap(other)
        for chromosome in self.genome_dic:
            if chromosome not in other.genome_dic:
                continue
            for frame in (-3, -2, -1, 1, 2, 3):
                self_it = self.genome_dic[chromosome][frame]
                other_it = other.genome_dic[chromosome][frame]
                overlap_it = IntervalTree()
                for interval in other_it.items():
                    for overlap_interval in self_it.overlap(interval):
                        overlap_it.add(overlap_interval)
                no_overlap_it = self_it.difference(overlap_it)
                no_overlap_ift.add_interval_tree(chromosome, frame, no_overlap_it)
        return no_overlap_ift

    def overlap_by_gene(self, other):
        """Get best intervals that overlaps a gene (set of exons)."""
        if not isinstance(other, IntervalFrameTree):
            raise TypeError("Input needs to be an interval frame tree.")

        if other.annotation_type != "gtf":
            raise TypeError("Only implemented for gtf file.")

        if set(self.genome_dic.keys()).intersection(set(other.genome_dic.keys())) == 0:
            raise KeyError("Both trees do not share any chromosomes.")

        gene_dic = {}
        for chromosome in self.genome_dic:
            if chromosome not in other.genome_dic:
                continue
            for frame in (-3, -2, -1, 1, 2, 3):
                self_it = self.genome_dic[chromosome][frame]
                other_it = other.genome_dic[chromosome][frame]
                for interval in other_it.items():
                    line_split = interval.data[1].split("\t")
                    info_line = " ".join(line_split[8:])
                    info_dic = info_line_to_dic(info_line)
                    gene_id = info_dic["gene_id"]
                    if gene_id not in gene_dic:
                        gene_dic[gene_id] = []
                    for overlap_interval in self_it.overlap(interval):
                        gene_dic[gene_id].append(overlap_interval.data[1])

        return gene_dic

    def make_annotation(self):
        """Build annoation from IFT.

        Only possible if IFT was build from Annotation and will return the same
        type of annotation as input.
        """
        if not self.annotation_type:
            raise ValueError("Interval annotation tree has no annotation type")
        lines = self.comment_lines
        for _chromosome, chromosome_dic in self.genome_dic.items():
            for _frame, interval_tree in chromosome_dic.items():
                lines += [i.data for i in interval_tree]
        lines.sort()
        return "".join([line[1] for line in lines])
