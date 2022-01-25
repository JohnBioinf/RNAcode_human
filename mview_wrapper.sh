#!/bin/bash

set -e

fasta=$1
html=$2

mview -in fasta -html data -bold -colormap CLUSTAL -coloring any "$fasta" -conservation on >> "$html"
