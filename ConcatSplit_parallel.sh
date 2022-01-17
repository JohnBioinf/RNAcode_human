#!/bin/bash

genome_alignment_dir=$1
min_size=$2
min_length=$3
max_len_no_split=$4
num_cpus=$5

out_file="$genome_alignment_dir/ConcatSplit_parallel.out"
log_file="$genome_alignment_dir/ConcatSplit_parallel.log"

echo "Start crunching" > "$out_file"
find "$genome_alignment_dir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | \
	parallel -j 20 --joblog "$log_file" \
    "python3 ConcatSplit.py $genome_alignment_dir/{} {} $min_size $min_length $max_len_no_split $num_cpus"

echo "Finished crunching" >> "$out_file"
