#!/bin/bash

genome_alignment_dir=$1
num_cpus=$2

out_file="$genome_alignment_dir/maf_stream_parallel.out"
log_file="$genome_alignment_dir/maf_stream_parallel.log"

echo "Start crunching" > "$out_file"
find "$genome_alignment_dir" -name "*.maf.gz" | \
	parallel -j "$num_cpus" --joblog "$log_file" \
    "python3 stream_chromosome.py {}"

echo "Finished crunching" >> "$out_file"
