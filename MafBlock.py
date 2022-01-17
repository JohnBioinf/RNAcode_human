#!/usr/bin/python3
"""Python module for class MafBlock."""


from copy import deepcopy
import subprocess


class MafBlock:
    """A maf block."""

    def __init__(self, lines=None):
        """Construct class."""
        self.block = []
        if lines:
            self.add(lines)

    def add(self, lines):
        """Add line or lines to maf-block."""
        if isinstance(lines, str):
            self.__add_string(lines)
        elif isinstance(lines, list):
            self.__add_list(lines)
        else:
            raise TypeError("Lines is not a string nor a list.")

    def __add_string(self, lines):
        if "\n" in lines:
            for entry in lines.split("\n"):
                if entry == "":
                    continue
                if entry[0] not in ["s", "a"]:
                    continue
                self.block.append(entry.split())
        else:
            if lines[0] in ["s", "a"]:
                self.block.append(lines.split())

    def __add_list(self, lines):
        self.block.append(lines)

    def __str__(self):
        """Show the maf-block formatted like a maf-file."""
        if self.is_empty():
            return "Empty maf\n\n"
        return "\n".join([" ".join(entry) for entry in self.block]) + "\n\n"

    def is_empty(self):
        """Test if object is empty."""
        return self.block == []

    def size(self):
        """Get the number of aligned sequences."""
        if self.is_empty():
            return 0
        return len(self.block)

    def __len__(self):
        """Get the length of the alignment."""
        if self.is_empty():
            return 0
        else:
            return len(self.block[1][-1])

    def coordinates(self):
        """Return start and end of maf-block."""
        if self.is_empty():
            return -1, -1
        return int(self.block[1][2]), int(self.block[1][2]) + int(self.block[1][3])

    def add_suffix(self, suffix, only_target=True):
        """Add suffix to names in maf block."""
        if self.is_empty():
            return
        if only_target:
            self.block[1][1] = self.block[1][1] + suffix
        else:
            for entry in self.block:
                if entry[0] == "a":
                    continue
                entry[1][1] = entry[1][1] + suffix

    def get_line_of_species(self, species):
        """Get index of line which is identified by species."""
        # return [i for i, line in enumerate(self.blocks) if line[1] == species][0]
        return [line for line in self.block if line[1] == species][0]

    def get_all_species(self, no_target=False):
        """Get all species in maf."""
        if self.is_empty():
            return []
        if no_target:
            start = 2
        else:
            start = 1
        return [line[1] for line in self.block[start:]]

    def _system_call(self, call_str, env=None, shell=False):
        """Perform system call based on input str.

        :param str call_str: The command as a string.
        :return: The return code and the output
        :rtype: tuple
        """
        with subprocess.Popen(
            call_str.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
            shell=shell,
        ) as process:
            process.wait()
            return process.returncode, process.communicate()[0].decode("UTF-8")

    def generate_html(self, out_dir, name=None, url_base=None, tmp_fasta_path="./tmp.fasta"):
        """Write out a nice alignment view of the maf block.

        Based on mview. (https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.67/mview-1.67.tar.gz)
        """
        if self.is_empty():
            print("Maf Block is empty.")
            return
        return_code, out = self._system_call("command -v mview")
        if return_code != 0:
            raise SystemError("mview is not available")

        with open(tmp_fasta_path, "w", encoding="UTF-8") as f_handle:
            for entry in self.block:
                if entry[0] == "a":
                    continue
                f_handle.write(f">{entry[1]}-{entry[2]}:{int(entry[3])+int(entry[2])}\n{entry[-1]}\n")
        if not name:
            name = self.block[1][1] if self.block[0][0] == "a" else self.block[0][1]

        call = f"./mview_wrapper.sh {tmp_fasta_path} {out_dir}/{name}.html"
        self._system_call(call)
        if url_base:
            print(f"See alignment under:\n{url_base}/{name}.html")

    def split_maf_block(self, split_start, split_end, suffix=""):
        """Split maf-block based on start and stop.

        The function will also change the genomic coordinates accordingly.

        :param list maf_block: list of list represtation of maf block.
        :param int split_start: start of split.
        :param int split_end: end  of split.
        :param str suffix: opional suffix that will be appended to name.
        :raise: indexError if split_start is smaller than 0 or if split_end is
            smaller or equal to split_start or if split_end is bigger than length
            of alignment.
        :return: smaller maf-block
        :rtype: list
        """
        if split_start < 0 or split_start > split_end:
            raise IndexError("Maf split start must be above 0.")
        new_maf = MafBlock()
        for entry in self.block:
            if entry[0] == "a":
                new_maf.add(entry)
                continue
            maf_type = entry[0]
            if suffix != "":
                name = entry[1] + "-" + str(suffix)
            else:
                name = entry[1]
            start = entry[2]
            size = entry[3]
            strand = entry[4]
            scr_size = entry[5]
            sequence = entry[6]
            if split_end > len(sequence):
                raise IndexError("Maf split end can not be bigger than sequence.")
            num_gaps_before_start = sequence[0:split_start].count("-")
            sequence = sequence[split_start:split_end]
            num_gaps = sequence.count("-")
            size = str(split_end - split_start - num_gaps)
            if strand == "+":
                start = str(int(start) + int(split_start) - num_gaps_before_start)
            else:
                start = str(int(start) + int(split_start) - num_gaps_before_start)
            new_maf.add(
                f"{maf_type} {name} {start} {size} {strand} {scr_size} {sequence}"
            )
        return new_maf

    def preprocess_block(self, max_len_no_split):
        """Preprocess bock by splitting long blocks into list of smaller blocks."""
        if len(self.block[1][-1]) < max_len_no_split:
            return [self]

        split_mafs = []
        end_reached = False
        for i, start in enumerate(
            range(0, len(self.block[1][-1]), int(max_len_no_split / 2))
        ):
            end = start + max_len_no_split
            if end >= len(self.block[1][-1]):
                end = len(self.block[1][-1])
                end_reached = True
            split_mafs += [self.split_maf_block(start, end, suffix=f"split-{i}")]
            if end_reached:
                break
        return split_mafs

    def add_symbols(self, symbol, positions):
        """Add symbol to every position in alignment."""
        for entry in self.block:
            if entry[0] == "a":
                continue
            for position in positions:
                entry[-1] = entry[-1][:position] + symbol + entry[-1][position:]

    def is_continuous_with(self, other):
        """Test if two maf blocks are continuous, only accept non existing species."""
        if self.is_empty():
            return False, -1

        if int(self.block[1][2]) > int(other.block[1][2]):
            raise ValueError("Maf block self is after maf block other")

        max_dist = 0
        discontious = False
        discontious_count = 0
        for species in set(other.get_all_species(no_target=True)).intersection(
            set(self.get_all_species(no_target=True))
        ):
            species_line_self = self.get_line_of_species(species)
            end_self = int(species_line_self[2]) + int(species_line_self[3])
            strand_self = species_line_self[4]

            species_line_other = other.get_line_of_species(species)
            start_other = int(species_line_other[2])
            strand_other = species_line_other[4]
            max_dist = max(start_other - end_self, max_dist)
            if (
                strand_self != strand_other
                or end_self + 12 < start_other
                or end_self > start_other
            ):
                # print(f"Dicontinous because of {species}-> strand:{strand_self != strand_other}, distance:{start_other - end_self}, reverse:{end_self > start_other}")
                discontious = True
                discontious_count += 1
        if discontious:
            print(f"Dicontinous because of {discontious_count}")
            return False, -1
        return True, max_dist

    def concat(self, other):
        """Concat to maf blocks."""
        continuous, max_dist = self.is_continuous_with(other)

        if not continuous:
            raise ValueError("Maf block self and other are not continuous")
        new_other = deepcopy(other)
        # Below could be a private method
        length_other = len(other)
        length_self = len(self)
        other_species = other.get_all_species(no_target=True)
        self_species = self.get_all_species(no_target=True)
        # if set(other_species) != set(self_species):
        #     print(f"Add species:{set(other_species).symmetric_difference(set(self_species))}")
        # if max_dist != 0:
        #     print("Add bridge")
        for species in self_species:
            if species in other_species:
                continue
            species_line_self = self.get_line_of_species(species)
            start = str(int(species_line_self[2]) + int(species_line_self[3]))
            strand = species_line_self[4]
            genome_length = species_line_self[5]
            new_other.add(
                ["s", species, start, "0", strand, genome_length, "-" * length_other]
            )

        for species in other_species:
            if species in self_species:
                continue
            species_line_other = other.get_line_of_species(species)
            genome_length = species_line_other[5]
            strand = species_line_other[4]
            start = species_line_other[2]
            self.add(
                ["s", species, start, "0", strand, genome_length, "-" * length_self]
            )

        for i, self_line in enumerate(self.block[1:]):
            if i == 0:
                other_line = new_other.block[1]
            else:
                species = self_line[1]
                other_line = new_other.get_line_of_species(species)
            end_self = int(self_line[2]) + int(self_line[3])
            start_other = int(other_line[2])
            dist = start_other - end_self
            gap_seq = dist * "X" + (max_dist - dist) * "-"
            # new length
            self_line[3] = str(int(self_line[3]) + int(other_line[3]) + dist)
            self_line[6] = self_line[6] + gap_seq + other_line[6]

#    def concat(self, other):
#        """Concat to maf blocks."""
#        if self.block == []:
#            self.block = other.block
#            return
#        if other.block == []:
#            return
#        new_other = deepcopy(other)
#        other_species = other.get_all_species()
#        self_species = self.get_all_species()
#        for species in self_species:
#            if species in other_species:
#                continue
#            overhang_nucs = self.get_overhanging_nucs(species)
#            genome_length = [line[5] for line in other.block if line[1] == species][0]
#            seq = "X" * overhang_nucs + (len(other) - overhang_nucs) * "-"
#            new_other.add(["s", species, 0, 0, "+", genome_length, seq])
#
#        for species in other_species:
#            if species in other_species:
#                continue
#            genome_length = [line[5] for line in other.block if line[1] == species][0]
#            self.add(["s", species, 0, 0, "+", genome_length, len(other) * "-"])
#
#        for species in set(self_species, other_species):
#            self_species_index = self.get_line_of_species(species)
#            other_species_index = new_other.get_line_of_species(species)
#            self.block[self_species_index][-1] = self.block[self_species_index][-1] + new_other.block[other_species_index][-1]
#
#    def get_overhanging_nucs(self, species):
#        """Get the number of nucleotides that should overhang into next block to preserve frame."""
#        if self.block[0][0] == "a":
#            target_seq = self.block[1][-1]
#        else:
#            target_seq = self.block[0][-1]
#        species_seq = [line[5] for line in self.block if line[1] == species][0]
#        nuc_count_target = 0
#        nuc_count_species = 0
#        for i in range(len(target_seq) - 1, -1, -1):
#            if target_seq[i] == "-":
#                target_gap = True
#            else:
#                target_gap = False
#                nuc_count_target += 1
#            if species_seq[i] == "-":
#                species_gap = True
#            else:
#                species_gap = False
#                nuc_count_species += 1
#            if not target_gap and not species_gap:
#                break
#        else:
#            raise ValueError(
#                f"No nucleotide alignes between target and {species} in {self.block[1][1]}"
#            )
#        nuc_diff = abs(nuc_count_species - nuc_count_species)
#        return nuc_diff % 3
