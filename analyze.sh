#!/bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_rt=12:00:00
#$ -l virtual_free=1
#$ -l ironfs
#$ -l gpu

source ~cs86/students/bin/PyRosetta.ScientificLinux-r53335.64Bit/SetPyRosettaEnvironment.sh
python analyze.py 1le0_m1.pdb output1.pdb 
