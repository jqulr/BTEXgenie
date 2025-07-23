#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

#
# diamond.sh — combine archetype sequences along with homogolous sequences identified through DIAMOND search
#
# Description:
#   For each input cluster FASTA, this script runs DIAMOND to pull in homologs
#   from NR, merges them with the original cluster, and writes out a combined
#   “raw” FASTA. 
#

# Load conda and environment
source /home/carmbrus/miniconda3/etc/profile.d/conda.sh
conda activate /home/juneq/.conda/envs/hmm

# Paths
NR_FASTA="/home/juneq/databases/nr.faa"
HG_DIR="/home/juneq/hmm/clusters_all_seqs_cleaned-split_fastas"
OUTDIR="/home/juneq/hmm/DIAMOND_HITS"
THREADS=16

mkdir -p "$OUTDIR"

for hg in "$HG_DIR"/*.faa; do
    # strip the .faa → "1-seqs_in_cluster-AKU14310_1"
    tag=$(basename "$hg" .faa)

    # if we've already created the merged _raw.faa, skip
    if [[ -s "$OUTDIR/${tag}_raw.faa" ]]; then
        echo "[SKIP] ${tag} already processed."
        continue
    fi

    echo "[INFO] Processing cluster $tag"

    # 1) DIAMOND search
    diamond blastp \
        --db "${NR_FASTA%.faa}.dmnd" \
        --query "$hg" \
        --evalue 1e-4 \
        --query-cover 70 \
        --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore \
        --out "$OUTDIR/${tag}.m8" \
        --threads "$THREADS" \
        --max-target-seqs 0

    # 2) Extract unique subject IDs
    cut -f2 "$OUTDIR/${tag}.m8" | sort -u > "$OUTDIR/${tag}_acc.txt"

    # 3) Fetch sequences if any
    if [[ -s "$OUTDIR/${tag}_acc.txt" ]]; then
        echo "[INFO] Fetching $(wc -l < "$OUTDIR/${tag}_acc.txt") sequences for $tag"
        seqtk subseq \
            "$NR_FASTA" \
            "$OUTDIR/${tag}_acc.txt" \
          > "$OUTDIR/${tag}_hits.faa"
    else
        echo "[WARN] No hits for $tag; touching empty file"
        : > "$OUTDIR/${tag}_hits.faa"
    fi

    # 4) Merge original cluster + hits
    cat "$hg" "$OUTDIR/${tag}_hits.faa" \
      > "$OUTDIR/${tag}_raw.faa"
done

echo "[INFO] All clusters processed."
