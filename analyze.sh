#!/bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_rt=12:00:00
#$ -l virtual_free=1
#$ -l ironfs
#$ -l gpu

source ~cs86/students/bin/PyRosetta.ScientificLinux-r53335.64Bit/SetPyRosettaEnvironment.sh
python analyze.py 4xdx.clean.pdb output_4xdx.clean_centroid_10000.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_centroid_1000.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_centroid_100.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_centroid_10.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_fullatom_10000.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_fullatom_1000.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_fullatom_100.pdb
python analyze.py 4xdx.clean.pdb output_4xdx.clean_fullatom_10.pdb

python analyze.py 4ytv.clean.pdb output_4ytv.clean_centroid_10000.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_centroid_1000.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_centroid_100.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_centroid_10.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_fullatom_10000.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_fullatom_1000.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_fullatom_100.pdb
python analyze.py 4ytv.clean.pdb output_4ytv.clean_fullatom_10.pdb
