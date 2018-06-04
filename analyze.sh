#!/bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_rt=1:00:00
#$ -l virtual_free=1
#$ -l ironfs
#$ -l gpu

source ~cs86/students/bin/PyRosetta.ScientificLinux-r53335.64Bit/SetPyRosettaEnvironment.sh
python analyze.py 4yud.clean.pdb output_4yud.clean_centroid_10000.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_centroid_1000.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_centroid_100.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_centroid_10.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_fullatom_10000.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_fullatom_1000.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_fullatom_100.pdb
python analyze.py 4yud.clean.pdb output_4yud.clean_fullatom_10.pdb

python analyze.py 5w0h.clean.pdb output_5w0h.clean_centroid_10000.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_centroid_1000.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_centroid_100.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_centroid_10.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_fullatom_10000.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_fullatom_1000.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_fullatom_100.pdb
python analyze.py 5w0h.clean.pdb output_5w0h.clean_fullatom_10.pdb

python analyze.py 6elm.clean.pdb output_6elm.clean_centroid_10000.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_centroid_1000.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_centroid_100.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_centroid_10.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_fullatom_10000.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_fullatom_1000.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_fullatom_100.pdb
python analyze.py 6elm.clean.pdb output_6elm.clean_fullatom_10.pdb

python analyze.py 1le0.clean.pdb output_1le0.clean_centroid_10000.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_centroid_1000.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_centroid_100.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_centroid_10.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_fullatom_10000.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_fullatom_1000.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_fullatom_100.pdb
python analyze.py 1le0.clean.pdb output_1le0.clean_fullatom_10.pdb

python analyze.py 1ubq.clean.pdb output_1ubq.clean_centroid_10000.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_centroid_1000.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_centroid_100.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_centroid_10.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_fullatom_10000.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_fullatom_1000.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_fullatom_100.pdb
python analyze.py 1ubq.clean.pdb output_1ubq.clean_fullatom_10.pdb
