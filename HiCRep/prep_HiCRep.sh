#!/bin/bash
#$ -l arch=linux-x64
#$ -l h_rt=4:0:0
#$ -l mem_free=16G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-11

i=$(($SGE_TASK_ID - 1))
input_names=("IJ053" "MS081" \
"IJ052" "MS051" "MS082" \
"IJ050" "MS052" "MS083" \
"IJ051" "MS053" "MS084")

INPUT=${input_names[i]}

echo $INPUT

python /netapp/home/bingkun/scripts/reg2hicrep_new.py 105000 58585000 5000 chr19 ${INPUT}.chr19.merged.cut final_matrix/${INPUT}.chr19.matrix
