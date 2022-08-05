#!/bin/bash

# This is script is a wrapper script which calculates the recall rate for
# RNAcode compared to the standart annotation.

set -e
set -u

DATA_DIR="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['DATA_DIR'])" \
    < "../parameters_local.json")"

DATA_DIR_OLD="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['DATA_DIR_OLD'])" \
    < "../parameters_local.json")"

HG38_ANNOTATION_WEB_FTP="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['HG38_ANNOTATION_WEB_FTP'])" \
	< "../parameters_local.json")"

R_PATH="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['R_PATH'])" \
	< "../parameters_local.json")"

TMP_DIR="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['TMP_DIR'])" \
	< "../parameters_local.json")"

export R_LIBS_USER="$R_PATH/libs"
export PATH="$R_PATH/R/bin/:$PATH"
# shellcheck disable=SC1091
source ../venv_k61/bin/activate

genome_alignment_dir="$DATA_DIR/multiz100way/"
# genome_alignment_dir="$DATA_DIR_OLD/multiz100way/"

echo "Download and count exons and genes."
# Download and modify standard annotation. 
ANNOTATION_FILE_GZ="${HG38_ANNOTATION_WEB_FTP##*/}"
ANNOTATION_FILE="${ANNOTATION_FILE_GZ//.gz/}"
ANNOTATION_FILE_MOD="${ANNOTATION_FILE//.gtf/.mod.gtf}"
# wget -nv -P "$DATA_DIR" "$HG38_ANNOTATION_WEB_FTP"
# gzip -d "$DATA_DIR/$ANNOTATION_FILE_GZ"

# Modify chromsom headers (add "chr") and discard all non coding genes. 
# The script print the number of genes and exons. The script folows the same
# naming convention for the annotation files as in this script and according
# with the parameters_local.json
# python3 ./NumberExonsGenes.py

echo "Calculate overlap and complement"
frames=(-3 -2 -1 1 2 3)
printf "" > "$genome_alignment_dir/RNAcode_no_overlap.bed"
printf "" > "$genome_alignment_dir/RNAcode_overlap.bed"
for frame in "${frames[@]}"; do
    echo "$frame"
    # Take from both annotations files only the genes which are the given frame.
    ./SubFrame.py "$frame" "$DATA_DIR/$ANNOTATION_FILE_MOD" > "$TMP_DIR/gtf"
    ./SubFrame.py "$frame" "$genome_alignment_dir/RNAcode.bed" > "$TMP_DIR/bed"

    # Get intersect. True positives 
    bedtools intersect -s -u -a "$TMP_DIR/bed" -b "$TMP_DIR/gtf" > /tmp/bed_overlap
    cat /tmp/bed_overlap >> "$genome_alignment_dir/RNAcode_overlap.bed"

    # Get relative complement. False positives 
    bedtools intersect -s -v -a "$TMP_DIR/bed" -b "$TMP_DIR/gtf" > /tmp/bed_no_overlap
    cat /tmp/bed_no_overlap >> "$genome_alignment_dir/RNAcode_no_overlap.bed"
done

echo "Sum the recalls"
python3 ./SumRecall.py "$genome_alignment_dir"
echo "Calculate recall by gene (if at least one exon is predicted whole gene is called)"
python3 ./recall.py "$genome_alignment_dir"
echo "Produce plot"
Rscript ./fdr_recal.R "$genome_alignment_dir"
