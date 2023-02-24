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

echo "Download and count exons and genes."
# Download and modify standard annotation. 
ANNOTATION_FILE_GZ="${HG38_ANNOTATION_WEB_FTP##*/}"
ANNOTATION_FILE="${ANNOTATION_FILE_GZ//.gz/}"
ANNOTATION_FILE_PSEUDO="${ANNOTATION_FILE//.gtf/.pseudo.gtf}"
ANNOTATION_FILE_CDS="${ANNOTATION_FILE//.gtf/.cds.gtf}"
# wget -nv -P "$DATA_DIR" "$HG38_ANNOTATION_WEB_FTP"
gzip -d "$DATA_DIR/$ANNOTATION_FILE_GZ"

# Modify chromsom headers (add "chr") and discard all non coding genes. 
# The script print the number of genes and exons. The script folows the same
# naming convention for the annotation files as in this script and according
# with the parameters_local.json
python3 ./NumberExonsGenes.py

for data_dir in $DATA_DIR_OLD $DATA_DIR; do
    genome_alignment_dir="$data_dir/multiz100way/"
    echo "Calculate overlap and complement in $data_dir"

    bedtools intersect -s -v -a "$DATA_DIR/$ANNOTATION_FILE_CDS" -b "$genome_alignment_dir/RNAcode_clean_tab.bed" > "$genome_alignment_dir/all_found_exons.gtf"
    frames=(-3 -2 -1 1 2 3)
    printf "" > "$genome_alignment_dir/RNAcode_no_overlap.bed"
    printf "" > "$genome_alignment_dir/RNAcode_overlap.bed"
    for frame in "${frames[@]}"; do
        echo "$frame"
        # Take from both annotations files only the genes which are the given frame.
        python3 ./SubFrame.py "$frame" "$DATA_DIR/$ANNOTATION_FILE_CDS" > "$TMP_DIR/gtf_cds"
        python3 ./SubFrame.py "$frame" "$genome_alignment_dir/RNAcode_clean.bed" > "$TMP_DIR/bed"
        cat "$TMP_DIR/gtf_cds" "$DATA_DIR/$ANNOTATION_FILE_PSEUDO" > "$TMP_DIR/gtf_pseudo"

        # Get intersect. True positives 
        # Find those that overlap with any CDS gene
        bedtools intersect -s -u -a "$TMP_DIR/bed" -b "$TMP_DIR/gtf_cds" > /tmp/bed_overlap
        cat /tmp/bed_overlap >> "$genome_alignment_dir/RNAcode_overlap.bed"

        # Get relative complement. False positives 
        # Here pseudo genes are included as are probably false but might be detected
        # so they are neither a true nor a false positiv 
        bedtools intersect -s -v -a "$TMP_DIR/bed" -b "$TMP_DIR/gtf_pseudo" > /tmp/bed_no_overlap
        cat /tmp/bed_no_overlap >> "$genome_alignment_dir/RNAcode_no_overlap.bed"
    done

    # echo "Sum the recalls"
    python3 ./SumRecall.py "$genome_alignment_dir"
    echo "Calculate recall by gene (if at least one exon is predicted whole gene is called)"
    python3 ./recall.py "$genome_alignment_dir"
done

echo "Produce plot"
Rscript ./fdr_recal_compare.R
