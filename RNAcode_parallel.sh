#!/bin/bash

genome_alignment_dir=$1
grouping=$2
num_cpus=$3

if [[ "$grouping" == "big_blocks" ]]; then
    reg_ex='big_block_[0-9]*\.maf'
elif [[ "$grouping" == "single_blocks" ]]; then
    reg_ex='big_block_[0-9]*-s_[0-9]*\.maf'
else
    echo "Error grouping: $grouping undefined."
fi

out_file="$genome_alignment_dir/RNAcode_${grouping}_parallel.out"
log_file="$genome_alignment_dir/RNAcode_${grouping}_parallel.log"

echo "Start crunching" > "$out_file"
find "$genome_alignment_dir" -name "$reg_ex" | \
	parallel -j "$num_cpus" --joblog "$log_file" \
	"RNAcode -t -o {}.res.tsv {} 2> {}.err > {}.out" \
	&>> "$out_file"

num_failed_jobs=$(awk 'NR>1 {if($7 != 0){print $0}}' "$log_file" | wc -l)

if [[ $num_failed_jobs -gt 0 ]]; then
	echo "$num_failed_jobs are failed retry." >> "$out_file"
	parallel -j "$num_cpus" --joblog "$log_file" --retry-failed
fi

echo "Finished crunching" >> "$out_file"
