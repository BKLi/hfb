#!/bin/bash
#$ -l arch=linux-x64
#$ -l h_rt=1:0:0
#$ -l mem_free=6G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-11

export PATH=/netapp/home/bingkun/seqtk:$PATH

i=$(($SGE_TASK_ID - 1))
input_names=("IJ053" "MS081" \
"IJ052" "MS051" "MS082" \
"IJ050" "MS052" "MS083" \
"IJ051" "MS053" "MS084")

INPUT=${input_names[i]}

echo $INPUT

sed 1d reg_raw.chr19.${INPUT}_merged_trimmed100_fastp.5k.and.MAPS2_pospoisson > reg_raw.chr19.${INPUT}_merged_trimmed100_fastp.5k.and.MAPS2_pospoisson_noheader
sed 1d reg_raw.chr19.${INPUT}_merged_trimmed100_fastp.5k.xor.MAPS2_pospoisson > reg_raw.chr19.${INPUT}_merged_trimmed100_fastp.5k.xor.MAPS2_pospoisson_noheader

cat reg_raw.chr19.${INPUT}_merged_trimmed100_fastp.5k.and.MAPS2_pospoisson_noheader reg_raw.chr19.${INPUT}_merged_trimmed100_fastp.5k.xor.MAPS2_pospoisson_noheader > ${INPUT}.chr19.merged
