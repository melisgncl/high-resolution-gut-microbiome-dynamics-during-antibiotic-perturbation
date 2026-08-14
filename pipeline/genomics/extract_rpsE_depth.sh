#!/usr/bin/env bash
# Dump per-position read depth across the rpsE region (+/-200 bp) for M5/M6/M8 vs I6.
set -uo pipefail
cd /mnt/c/Users/melis/Desktop/bioinformatics-portfolio/hgt-study
if [ -f ~/miniforge3/etc/profile.d/conda.sh ]; then source ~/miniforge3/etc/profile.d/conda.sh
elif [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then source ~/miniconda3/etc/profile.d/conda.sh
fi
conda activate paeni-genomics
OUT=results/genomics/comparative
REGION="contig_1:228119-229016"   # rpsE 228319-228816 +/-200
for name in M5 M6 M8; do
  samtools depth -a -r "$REGION" "$OUT/${name}_reads_vs_I6.bam" > "$OUT/${name}_rpsE_depth.tsv"
done
echo "DONE_DEPTH"
