#!/bin/bash
set -euo pipefail

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate paeni-genomics

PROJ="/mnt/c/Users/melis/Desktop/bioinformatics-portfolio/hgt-study"
SEQD="$PROJ/data/sequence_data/TPHP3P_results"
REF="$PROJ/data/references/GCF_022494515.1_genomic.fna"
WORK="$HOME/fastani_work"
OUTD="$PROJ/results/genomics/comparative"

mkdir -p "$WORK"

# Build query list (M5, M6, M8 — skip M7)
echo "$SEQD/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna" >  "$WORK/query_list.txt"
echo "$SEQD/TPHP3P_2_M6__+_/annotation/TPHP3P_2_M6__+_.fna" >> "$WORK/query_list.txt"
echo "$SEQD/TPHP3P_4_M8__+_/annotation/TPHP3P_4_M8__+_.fna" >> "$WORK/query_list.txt"

echo "=== FastANI query list ==="
cat "$WORK/query_list.txt"

echo ""
echo "=== Running FastANI ==="
fastANI --ql "$WORK/query_list.txt" -r "$REF" -o "$WORK/fastani_out.txt" 2>"$WORK/fastani.err"
echo "FastANI exit: $?"

echo ""
echo "=== FastANI results ==="
cat "$WORK/fastani_out.txt"

# Copy results to project output dir
cp "$WORK/fastani_out.txt" "$OUTD/fastani_M5M6M8_vs_I6.txt"
echo ""
echo "Copied to $OUTD/fastani_M5M6M8_vs_I6.txt"
