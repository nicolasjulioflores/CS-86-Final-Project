#!/bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_rt=12:00:00
#$ -l virtual_free=1
#$ -l ironfs
#$ -l gpu

source ~cs86/students/bin/setupPyr.sh
python struc_predict_repack.py