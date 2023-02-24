#!/usr/bin/python3
"""Python module for class MafBlock."""


from copy import deepcopy
import subprocess
import gzip
import re
from collections import defaultdict


def append_html(string, name):
    """Do maf block."""
    with open(
        f"/homes/biertruck/john/public_html/mview/{name}.html", "a", encoding="UTF-8"
    ) as f_handle:
        f_handle.write(string)


def percent_undefined(string):
    """Calculate the percentage of gaps and undefined nucs (X)."""
    return (len(re.findall("[Nn-]", string)) / len(string)) * 100


class MafBlock:
    """A maf block."""

    def __init__(self, lines=None, allowed_dist=12):
        """Construct class."""
        self.block = []
        self.allowed_dist = allowed_dist
        self.block_index_list = []
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

    def add_index(self, index_object):
        """Add index to block index list."""
        if isinstance(index_object, list):
            for index in index_object:
                self.add_index(index)
        elif isinstance(index_object, int):
            self.block_index_list.append(index_object)
        else:
            raise ValueError("add_index() needs list or int.")

    def sort(self):
        """Sort maf block based on species name."""
        self.block[:1] = sorted(self.block[:1], key=lambda x: x[1])

    def is_empty(self):
        """Test if object is empty."""
        return self.block == []

    def size(self):
        """Get the number of aligned sequences."""
        if self.is_empty():
            return 0
        return len(self.block) - 1

    def __len__(self):
        """Get the length of target sequence with gaps."""
        if self.is_empty():
            return 0
        else:
            return len(self.block[1][-1])

    def len_no_gaps(self):
        """Get the length of target sequence without gaps."""
        if self.is_empty():
            return 0
        else:
            return len(self.block[1][-1].replace("-", ""))

    def get_target(self):
        """Get name of target species."""
        return self.block[1][1]

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

    def delete_species(self, species):
        """Delte species from block."""
        if isinstance(species, list):
            for species_name in species:
                self.delete_species(species_name)
        else:
            self.block = [line for line in self.block if line[1] != species]

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

    def generate_html(
        self, out_dir, name=None, url_base=None, tmp_fasta_path="./tmp.fasta"
    ):
        """Write out a nice alignment view of the maf block.

        Based on mview. (https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.67/mview-1.67.tar.gz)
        """
        html_file_path = f"{out_dir}/{name}.html"
        if self.is_empty():
            with open(html_file_path, "a", encoding="UTF-8") as f_handle:
                f_handle.write("<p>Empty Maf Block</p><br>")
            return
        self.sort()
        return_code, out = self._system_call("command -v mview")
        if return_code != 0:
            raise SystemError("mview is not available")

        with open(tmp_fasta_path, "w", encoding="UTF-8") as f_handle:
            for entry in self.block:
                if entry[0] == "a":
                    continue
                f_handle.write(
                    f">{entry[1]}_{entry[2]}:{int(entry[3])+int(entry[2])}\n{entry[-1]}\n"
                )
        if not name:
            name = self.block[1][1] if self.block[0][0] == "a" else self.block[0][1]

        call = f"./mview_wrapper.sh {tmp_fasta_path} {html_file_path}"
        return_code, out = self._system_call(call)
        if return_code != 0 or out != "":
            print(out)
            raise SystemError("mview exited with error")

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
        """Test if two maf blocks are continuous.

        Coninuity is based on the distance and the orientation of the sequences
        in the genome. If the distance between two sequences of the same species
        in two maf blocks is greater than the attribute allowed_dist the two
        blocks are discontinous. If the strand changes or the next block starts
        before the current block the two blocks are also discontinous. Ignored
        are species that only exist in one block but not the other, because
        these species can be easily replaced with only gaps in the next block.

        :param MafBlock arg1: Another MafBlock.
        :return: bigest distance of two sequences between any species which are
            continous. And a list of species which are discontinous.
        :rtype: (int, list)
        """
        if int(self.block[1][2]) > int(other.block[1][2]):
            raise ValueError("Maf block self is after maf block other")

        other_species = set(other.get_all_species())
        self_species = set(self.get_all_species())

        max_dist = 0
        discontious_species = []
        for species in other_species.intersection(self_species):
            species_line_self = self.get_line_of_species(species)
            end_self = int(species_line_self[2]) + int(species_line_self[3])
            strand_self = species_line_self[4]

            species_line_other = other.get_line_of_species(species)
            start_other = int(species_line_other[2])
            strand_other = species_line_other[4]
            # print(species)
            # print(start_other - end_self)
            if (
                strand_self != strand_other
                or end_self + self.allowed_dist < start_other
                or end_self > start_other
            ):
                # print(
                #     f"Discontious because of {species}-> strand:{strand_self != strand_other}, distance:{start_other - end_self}, reverse:{end_self > start_other}"
                # )
                discontious_species.append(species)
            else:
                max_dist = max(start_other - end_self, max_dist)

        return max_dist, discontious_species

    def concat(self, other, max_del=0):
        """Concat to maf blocks.

        :param MafBlock arg1: Another MafBlock.
        :param int arg1: The number of species that are allowed to be deleted.
        :return: 0 if concatination was sucessfull, negativ int if x species
            imped concatination
        :rtype: int
        """
        max_dist, discontious_species = self.is_continuous_with(other)
        # VERBOSE
        # print(f"The two blocks have a maximal distance {max_dist} and {discontious_species} discontious species.")

        if self.get_target() in discontious_species:
            return -float("inf")

        if len(discontious_species) > max_del:
            # VERBOSE
            # print(discontious_species)
            return -len(discontious_species)

        new_other = deepcopy(other)

        self.add_index(new_other.block_index_list)

        # This deletes species if no species imped, list should be empty.
        self.delete_species(discontious_species)
        new_other.delete_species(discontious_species)

        # Below could be a private method
        length_other = len(other)
        length_self = len(self)
        other_species = set(new_other.get_all_species(no_target=True))
        self_species = set(self.get_all_species(no_target=True))

        # VERBOSE
        # if other_species != self_species:
        #     print(f"Add species to other: {self_species - other_species}")
        #     print(f"Add species to self: {other_species - self_species}")
        #     # print(other_species)
        #     # print(self_species)
        # if max_dist != 0:
        #     print(f"Add bridge length {max_dist}")

        # Add new species to other
        for species in self_species - other_species:
            species_line_self = self.get_line_of_species(species)
            start = str(int(species_line_self[2]) + int(species_line_self[3]))
            strand = species_line_self[4]
            genome_length = species_line_self[5]
            new_other.add(
                ["s", species, start, "0", strand, genome_length, "-" * length_other]
            )

        # Add new species to self
        for species in other_species - self_species:
            species_line_other = other.get_line_of_species(species)
            genome_length = species_line_other[5]
            strand = species_line_other[4]
            start = species_line_other[2]
            self.add(
                ["s", species, start, "0", strand, genome_length, "-" * length_self]
            )

        # Concat
        for i, self_line in enumerate(self.block[1:]):
            # First line is target
            if i == 0:
                other_line = new_other.block[1]
            # Get the other line
            else:
                species = self_line[1]
                other_line = new_other.get_line_of_species(species)
            end_self = int(self_line[2]) + int(self_line[3])
            start_other = int(other_line[2])
            dist = start_other - end_self
            gap_seq = dist * "N" + (max_dist - dist) * "-"
            # new length
            self_line[3] = str(int(self_line[3]) + int(other_line[3]) + dist)
            # new seq
            self_line[6] = self_line[6] + gap_seq + other_line[6]

        return 0

    def clean(self):
        """Clean maf block.

        Cut rows with only gaps in target at beginning and end.
        Remove rows which only consists of rows.
        """
        try:
            start_row = [
                m.start() + len(m.group()) for m in re.finditer("^-+", self.block[1][-1])
            ][0]
        except IndexError:
            start_row = 0

        try:
            end_row = [m.start() for m in re.finditer("-+$", self.block[1][-1])][0]
        except IndexError:
            end_row = len(self.block[1][-1])

        self.block[1:] = [
            line[:-1] + [line[-1][start_row:end_row]] for line in self.block[1:]
        ]

        self.block[2:] = [
            line for line in self.block[2:] if percent_undefined(line[-1]) < 95
        ]


class MafStream:
    """A generator that iterates over a maf file."""

    def __init__(
        self,
        path=None,
        min_length_del=None,
        max_del_species=None,
        min_size=None,
        min_length=None,
        max_len_no_split=None,
    ):
        """Init stream with path to file."""
        for arg, value in locals().items():
            if arg == "self":
                continue
            if value is None:
                raise ValueError("all arguments must be set")

            setattr(self, arg, value)
        # Saves off set for each block to make iteration faster.
        # key: block_index val: off_set
        self.off_set_dic = {}
        # Previous calculated blocks
        # key: block_index val: MafBlock()
        self.block_dic = {}

    def __iter__(self):
        """Iterat over maf blocks in file."""
        if self.path.endswith("gz"):
            open_fn = gzip.open
            read_arg = "rb"
        else:
            open_fn = open
            read_arg = "r"
        block_index = 0
        maf = MafBlock()
        with open_fn(self.path, read_arg) as maf_handle:
            while True:
                line = maf_handle.readline()
                try:
                    line = line.decode("UTF8")
                except AttributeError:
                    pass
                # EOF
                if line == "":
                    break

                if line != "\n":
                    maf.add(line)
                else:
                    maf.add_index(block_index)
                    yield maf
                    maf = MafBlock()
                    block_index += 1

    def iterate_from_use_offset(self, block_index):
        """Iterate over maf blocks starting with specific block index.

        The attribute off_set_dic is used to make jumps.
        """
        if self.path.endswith("gz"):
            open_fn = gzip.open
            read_arg = "rb"
        else:
            open_fn = open
            read_arg = "r"

        if block_index in self.off_set_dic:
            current_block_index = block_index
            off_set = self.off_set_dic[current_block_index]
        # If block index is not in dic get highest block index
        # This should be used most of the time. Maybe a list would be quicker
        # compared to dic?
        elif len(self.off_set_dic) == 0:
            current_block_index = 0
            off_set = 0
        else:
            current_block_index = sorted(list(self.off_set_dic.keys()))[-1]
            off_set = self.off_set_dic[current_block_index]

        maf = MafBlock()
        with open_fn(self.path, read_arg) as maf_handle:
            maf_handle.seek(off_set)
            while True:
                line = maf_handle.readline()
                try:
                    line = line.decode("UTF8")
                except AttributeError:
                    pass
                # EOF
                if line == "":
                    if current_block_index < block_index:
                        raise IndexError("block index out of range")
                    break

                if line[0] == "a":
                    self.off_set_dic[current_block_index] = off_set

                off_set += len(line)

                if line != "\n":
                    maf.add(line)
                else:
                    maf.add_index(current_block_index)
                    if current_block_index >= block_index:
                        yield maf
                    maf = MafBlock()
                    current_block_index += 1

    def iterate_from(self, block_index):
        """Iterate over maf blocks starting with specific block index."""
        if self.path.endswith("gz"):
            open_fn = gzip.open
            read_arg = "rb"
        else:
            open_fn = open
            read_arg = "r"

        current_block_index = 0
        maf = MafBlock()
        with open_fn(self.path, read_arg) as maf_handle:
            while True:
                line = maf_handle.readline()
                try:
                    line = line.decode("UTF8")
                except AttributeError:
                    pass
                # EOF
                if line == "":
                    break

                if line != "\n":
                    maf.add(line)
                else:
                    maf.add_index(current_block_index)
                    if current_block_index >= block_index:
                        yield maf
                    maf = MafBlock()
                    current_block_index += 1

        if current_block_index < block_index:
            raise IndexError("block index out of range")

    def concat_blocks_trivial(self, position=(), index_range=(0, float("inf"))):
        """Concatinate blocks without any deletion."""
        maf = MafBlock()

        for next_maf in self.iterate_from(index_range[0]):
            # Only first iteration
            if maf.is_empty():
                maf = next_maf
                continue
            if position != ():
                start_maf, end_maf = maf.coordinates()
                if end_maf < position[0]:
                    continue
                if start_maf > position[1]:
                    break
            if max(maf.block_index_list) > index_range[1]:
                yield maf
                break
            return_value = maf.concat(next_maf)
            if return_value != 0:
                yield maf
                maf = next_maf
        if max(maf.block_index_list) < index_range[1]:
            yield maf

    def concat_blocks_with_deletion(self, position=(), index_range=(0, float("inf"))):
        """Concatinate blocks with deletion."""
        maf = MafBlock()

        for next_maf in self.concat_blocks_trivial(position=position, index_range=index_range):
            # Only first iteration
            if maf.is_empty():
                maf = next_maf
                continue
            if (
                maf.len_no_gaps() < self.min_length_del
                and next_maf.len_no_gaps() < self.min_length_del
            ):
                smallest_size = min([maf.size(), next_maf.size()])
                impeding_species = abs(maf.concat(next_maf))
                ratio = impeding_species / smallest_size
                if ratio < 0.1:
                    maf.concat(next_maf, max_del=impeding_species)
                    continue

            yield maf
            maf = next_maf
        if maf != next_maf:
            yield maf

    def split_stream(self, position=(), index_range=(0, float("inf"))):
        """Split blocks with sliding window."""
        for maf in self.concat_blocks_with_deletion(position=position, index_range=index_range):
            if len(maf) < self.max_len_no_split:
                yield maf
                continue
            end_reached = False
            for i, start in enumerate(
                range(0, len(maf), int(self.max_len_no_split / 2))
            ):
                end = start + self.max_len_no_split
                if end >= len(maf):
                    end = len(maf)
                    end_reached = True
                yield maf.split_maf_block(start, end, suffix=f"split-{i}")
                if end_reached:
                    break

    def discard_stream(self, position=(), index_range=(0, float("inf"))):
        """Sort out small mafs and trim mafs."""
        for maf in self.split_stream(position=position, index_range=index_range):
            if maf.size() - 1 < self.min_size or maf.len_no_gaps() < self.min_length:
                continue
            maf.clean()
            if maf.size() - 1 < self.min_size or maf.len_no_gaps() < self.min_length:
                continue
            yield maf

    def concat_blocks(self, block_index, only_block=True, split=False):
        """Concatinate blocks starting with a specific block index."""
        maf_list = []
        maf = MafBlock()

        # might be set if look into the future becomes reality
        calc_block_index = -1

        for block_index, next_maf in self.iterate_from(block_index):
            next_maf.add_suffix(f"_block-index_{block_index}")

            # This can be triggered if new block starts and the current index
            # was already calculated
            if block_index in self.block_dic and maf.is_empty():
                calc_block_index, maf = self.block_dic[block_index]
                if only_block:
                    return calc_block_index, maf
                continue

            # Only first iteration
            if maf.is_empty():
                maf = next_maf
                continue

            if block_index < calc_block_index:
                continue

            return_value = maf.concat(next_maf)
            # Return value is zero if concatenation was successful
            if return_value == 0:
                pass
            # The returned negative int says how many species imped a
            # concatenation.
            elif return_value >= -self.max_del_species:
                # If the current block or the next block are below threshold
                # force concatenation even with deleting species.
                try:
                    # Test if the current maf block is to small
                    if maf.len_no_gaps() < self.min_length_del:
                        maf.concat(next_maf, max_del=self.max_del_species)
                    # Test if the next maf block is to small
                    elif next_maf.len_no_gaps() < self.min_length_del:
                        # Test if the next block can also not be extend to
                        # become big enough.
                        end_block_index, future_maf = self.concat_blocks(block_index)
                        # Safe future maf if force concatenation is triggered in
                        # does not need to be calculated again.
                        if future_maf.len_no_gaps() < self.min_length_del:
                            maf.concat(next_maf, max_del=self.max_del_species)
                        else:
                            # Safe future maf if force concatenation is not
                            # triggered. Does not need to be calculated again.
                            self.block_dic[block_index] = end_block_index, future_maf
                            return_value = -self.max_del_species - 1

                    # If both blocks are sufficiently big keep current block and
                    # proceed with concatenation with next block.
                    else:
                        return_value = -self.max_del_species - 1
                # Catch EOF
                except IndexError as e:
                    if str(e) == "block index out of range":
                        return_value = -self.max_del_species - 1
                    else:
                        raise

            # Here concatenation ends. Either to many impeding species or
            # impeding species is below self.max_del_species but both maf
            # blocks, the current and the next are big enough.
            if return_value < -self.max_del_species:
                if only_block:
                    return block_index, maf
                # Catches lower bound size
                if (
                    maf.size() - 1 < self.min_size
                    and len(maf.block[1][-1].replace("-", "")) < self.min_length
                ):
                    maf = next_maf
                    continue
                if split:
                    maf_list += maf.preprocess_block(self.max_len_no_split)
                else:
                    maf_list.append(maf)
                if block_index in self.block_dic:
                    calc_block_index, maf = self.block_dic[block_index]
                else:
                    maf = next_maf

        if only_block:
            return block_index, maf
        maf_list.append(maf)
        return block_index, maf_list

    def concat_blocks_old(self, block_index, only_block=True, split=False):
        """Concatinate blocks starting with a specific block index."""
        maf_list = []
        maf = MafBlock()
        imped_dic = defaultdict(lambda: 0)

        for block_index, next_maf in self.iterate_from(block_index):
            next_maf.add_suffix(f"_block-index_{block_index}")

            # Only first iteration
            if maf.is_empty():
                maf = next_maf
                continue

            return_value = maf.concat(next_maf)
            imped_dic[return_value] += 1
            # Return value is zero if concatenation was successful
            if return_value == 0:
                pass
            # The returned negative int says how many species imped a
            # concatenation.
            elif return_value >= -self.max_del_species:
                # If the current block or the next block are below threshold
                # force concatenation even with deleting species.
                try:
                    # Test if the current maf block is to small
                    if maf.len_no_gaps() < self.min_length_del:
                        maf.concat(next_maf, max_del=self.max_del_species)
                    # Test if the next maf block is to small
                    elif next_maf.len_no_gaps() < self.min_length_del:
                        # Test if the next block can also not be extend to
                        # become big enough.
                        future_maf = self.concat_blocks(block_index)
                        if future_maf.len_no_gaps() < self.min_length_del:
                            maf.concat(next_maf, max_del=self.max_del_species)
                        else:
                            return_value = -self.max_del_species - 1

                    # If both blocks are sufficiently big keep current block and
                    # proceed with concatenation with next block.
                    else:
                        return_value = -self.max_del_species - 1
                # Catch EOF
                except IndexError as e:
                    if str(e) == "block index out of range":
                        return_value = -self.max_del_species - 1
                    else:
                        raise

            # Here concatenation ends. Either to many impeding species or
            # impeding species is below self.max_del_species but both maf
            # blocks, the current and the next are big enough.
            if return_value < -self.max_del_species:
                if only_block:
                    return maf
                # Catches lower bound size
                if (
                    maf.size() - 1 < self.min_size
                    and len(maf.block[1][-1].replace("-", "")) < self.min_length
                ):
                    maf = next_maf
                    continue
                if split:
                    maf_list += maf.preprocess_block(self.max_len_no_split)
                else:
                    maf_list.append(maf)
                maf = next_maf

        if only_block:
            return maf

        with open("./bar_data_imped_chr10.csv", "w", encoding="UTF-8") as f_handle:
            f_handle.write("\n".join([f"{n},{s}" for n, s in imped_dic.items()]))

        maf_list.append(maf)
        return maf_list

    def concat_blocks_verbose(self, block_index, only_block=True, split=False):
        """Concatinate blocks starting with a specific block index, with verbose."""
        maf_list = []
        maf = MafBlock()

        # might be set if look into the future becomes reality
        calc_block_index = -1

        for block_index, next_maf in self.iterate_from(block_index):
            next_maf.add_suffix(f"_block-index_{block_index}")

            if block_index > 56:
                break
            append_html(f"<h2>Iteration: {block_index}</h2>\n")

            # This can be triggered if new block starts and the current index
            # was already calculated
            if block_index in self.block_dic and maf.is_empty():
                calc_block_index, maf = self.block_dic[block_index]
                if only_block:
                    append_html("<h3>Block was already calculated.</h3>\n")
                    return calc_block_index, maf
                append_html("<h3>Start Jump ...</h3>\n")
                append_html("<h3>Future Maf Block</h3>\n")
                maf.generate_html(
                    "/homes/biertruck/john/public_html/mview/", name="test_concat_parts"
                )
                continue

            if block_index < calc_block_index:
                append_html("<h3>... still jumping ...</h3>\n")
                continue

            print()
            print(f"Iteration: {block_index}")
            print(f"Length maf: {len(maf)}")
            for line in maf.block[1:]:
                print(f"{line[1]}: {len(line[-1])}")
            print(f"Length next maf: {len(next_maf)}")
            for line in next_maf.block[1:]:
                print(f"{line[1]}: {len(line[-1])}")
            append_html("<h3>Current Maf Block</h3>\n")
            maf.generate_html(
                "/homes/biertruck/john/public_html/mview/", name="test_concat_parts"
            )
            append_html("<h3>Next Maf Block</h3>\n")
            next_maf.generate_html(
                "/homes/biertruck/john/public_html/mview/", name="test_concat_parts"
            )

            # Only first iteration
            if maf.is_empty():
                maf = next_maf
                continue

            return_value = maf.concat(next_maf)
            # Return value is zero if concatenation was successful
            if return_value == 0:
                # VERBOSE
                append_html("<h4>Concatenation was successful</h4>\n")
                print("Concatenated")
                print(f"Length maf: {len(maf)}")
                pass
            # The returned negative int says how many species imped a
            # concatenation.
            elif return_value >= -self.max_del_species:
                append_html(
                    "<h4>Concatenation failed but only one species imped.</h4>\n"
                )
                print("Concatenation failed but only one species imped.")
                # input("Wait")
                # If the current block or the next block are below threshold
                # force concatenation even with deleting species.
                try:
                    # Test if the current maf block is to small

                    # VERBOSE
                    print(
                        f"Size of current maf is to small {maf.len_no_gaps() < self.min_length_del}"
                    )
                    print(
                        f"Size of next maf is to small {next_maf.len_no_gaps() < self.min_length_del}"
                    )
                    # print(f"Size of concatenated next maf is to small {self.concat_blocks(block_index + 1).len_no_gaps() < self.min_length_del}")

                    if maf.len_no_gaps() < self.min_length_del:
                        # VERBOSE
                        append_html(
                            "<h4>Current block is to small concat anyways.</h4>\n"
                        )
                        print("Current block is to small concat anyways.")

                        maf.concat(next_maf, max_del=self.max_del_species)
                    # Test if the next maf block is to small
                    elif next_maf.len_no_gaps() < self.min_length_del:
                        # VERBOSE
                        append_html(
                            "<h4>Next block is to small, look into future.</h4>\n"
                        )
                        print("Next block is to small, look into future.")
                        # input("Wait")
                        append_html('<div style="margin-left:70px">\n')

                        # Test if the next block can also not be extend to
                        # become big enough.
                        end_block_index, future_maf = self.concat_blocks_verbose(
                            block_index
                        )
                        # VERBOSE
                        append_html("</div>\n")
                        append_html("<h1>Back from the future</h1>\n")
                        append_html(f"<h2>Iteration: {block_index}</h2>\n")
                        append_html("<h3>Current Maf Block</h3>\n")
                        maf.generate_html(
                            "/homes/biertruck/john/public_html/mview/",
                            name="test_concat_parts",
                        )
                        append_html("<h3>Next Maf Block</h3>\n")
                        next_maf.generate_html(
                            "/homes/biertruck/john/public_html/mview/",
                            name="test_concat_parts",
                        )
                        append_html("<h3>Future block</h3>\n")
                        future_maf.generate_html(
                            "/homes/biertruck/john/public_html/mview/",
                            name="test_concat_parts",
                        )

                        print("Back from the future.")
                        if future_maf.len_no_gaps() < self.min_length_del:
                            print(
                                "Even in the future block is to small. Force concatenation."
                            )
                            append_html(
                                "<h4>Even in the future block is to small. Force concatenation.</h4>\n"
                            )
                            maf.concat(next_maf, max_del=self.max_del_species)
                        else:
                            # Safe future maf if force concatenation is not
                            # triggered. Does not need to be calculated again.
                            self.block_dic[block_index] = end_block_index, future_maf
                            return_value = -self.max_del_species - 1
                            append_html(
                                f"<h4>Keys block_dic = {list(self.block_dic.keys())}</h4>\n"
                            )
                            append_html("<h4>In the future block is big enough.</h4>\n")
                            print("In the future block is big enough.")
                        # VERBOSE
                        # input("Wait")

                    # If both blocks are sufficiently big keep current block and
                    # proceed with concatenation with next block.
                    else:
                        return_value = -self.max_del_species - 1
                # Catch EOF
                except IndexError as e:
                    if str(e) == "block index out of range":
                        return_value = -self.max_del_species - 1
                    else:
                        raise
            else:
                append_html("<h4>Concatenation failed to many species imped.</h4>\n")

            # Here concatenation ends. Either to many impeding species or
            # impeding species is below self.max_del_species but both maf
            # blocks, the current and the next are big enough.
            if return_value < -self.max_del_species:
                print("Concatination ends.")
                append_html("<h2>Concatination ends.</h4>\n")
                if only_block:
                    return block_index, maf
                # Catches lower bound size
                if (
                    maf.size() - 1 < self.min_size
                    and len(maf.block[1][-1].replace("-", "")) < self.min_length
                ):
                    maf = next_maf
                    continue
                if split:
                    maf_list += maf.preprocess_block(self.max_len_no_split)
                else:
                    maf_list.append(maf)
                print("Over write current maf")
                print(f"Length maf: {len(maf)}")
                print(f"Length next maf: {len(next_maf)}")
                if block_index in self.block_dic:
                    append_html("<h3>Start Jump ...</h3>\n")
                    calc_block_index, maf = self.block_dic[block_index]
                    append_html("<h3>Future Maf Block</h3>\n")
                    maf.generate_html(
                        "/homes/biertruck/john/public_html/mview/",
                        name="test_concat_parts",
                    )
                else:
                    maf = next_maf

        if only_block:
            return block_index, maf
        maf_list.append(maf)
        return block_index, maf_list
