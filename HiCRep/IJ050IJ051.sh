#!/bin/bash
#$ -l h_rt=20:0:0
#$ -l mem_free=60G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -m ae
#$ -M libingkun1997@gmail.com

Rscript IJ050IJ051.R IJ050.chr19.matrix IJ051.chr19.matrix