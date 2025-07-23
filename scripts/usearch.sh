#!/bin/bash
set -euo pipefail
IFS=$'\n\t'


#
# usearch.sh — Dereplicate and cluster DIAMOND hit FASTAs at 98% identity
#
# Description:
#   This script functions to remove duplicate sequences and cluster at 98% similarity
#

cd /home/juneq/hmm

input_dir="/home/juneq/hmm/DIAMOND_HITS"

output_dir="/home/juneq/hmm/DIAMOND_HITS/processed"
mkdir -p "$output_dir"

for file in "$input_dir"/*_raw.faa; do
    base=$(basename "$file" "_raw.faa")

    echo "processing $base"

    output="${output_dir}/${base}_uniq.fa"
    centroids="${output_dir}/${base}_centroids.fa"
    clusters_uc="${output_dir}/${base}.uc"

    # dereplicate duplicates
    /home/juneq/programs/usearch11.0.667_i86linux32  -fastx_uniques "$file" -fastaout "$output" -sizeout

    # cluster at 98% identity
    /home/juneq/programs/usearch11.0.667_i86linux32 -cluster_fast "$output" -id 0.98 -centroids "$centroids" -uc "$clusters_uc"

    echo "$base done processing."

done

echo "USEARCH completed."
