#!/usr/bin/env bash
# Map M5/M6/M8 Nanopore reads to the I6 reference and quantify the rpsE deletion.
# Run in WSL Ubuntu (conda env: paeni-genomics).
set -uo pipefail
cd /mnt/c/Users/melis/Desktop/bioinformatics-portfolio/hgt-study

if [ -f ~/miniforge3/etc/profile.d/conda.sh ]; then source ~/miniforge3/etc/profile.d/conda.sh
elif [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then source ~/miniconda3/etc/profile.d/conda.sh
fi
conda activate paeni-genomics

REF=data/references/bakta_I6/I6_bakta.fna       # I6 genome; rpsE = contig_1:228319-228816
OUT=results/genomics/comparative
RPSE="contig_1:228319-228816"
LOG=$OUT/rpsE_vs_I6_readlevel.log
: > "$LOG"
samtools faidx "$REF"

for S in M5:TPHP3P_1_M5 M6:TPHP3P_2_M6 M8:TPHP3P_4_M8; do
  name=${S%%:*}; pref=${S##*:}
  reads=data/sequence_data/TPHP3P_fastq/${pref}__+_.fastq.gz
  bam=$OUT/${name}_reads_vs_I6.bam
  echo "=== $name ===" | tee -a "$LOG"
  minimap2 -ax map-ont "$REF" "$reads" 2>>"$LOG" | samtools sort -o "$bam" -
  samtools index "$bam"
  samtools depth -a -r "$RPSE" "$bam" \
    | awk '{s+=$3; if(NR==1||$3<mn)mn=$3; if($3>mx)mx=$3; n++} END{printf "  cov mean=%.1f min=%d max=%d\n",s/n,mn,mx}' | tee -a "$LOG"
  total=$(samtools view -F 0x900 "$bam" "$RPSE" | wc -l)
  del9=$(samtools view -F 0x900 "$bam" "$RPSE" | awk '$6 ~ /9D/' | wc -l)
  echo "  primary reads over rpsE: $total ; carrying 9D deletion: $del9" | tee -a "$LOG"
done
echo "DONE_RPSE_VS_I6"
