from rosetta import *
import sys

init()

'''
Usage:
analyze.py clean_pdb repacked_pdb
'''

argc = len(sys.argv)

if argc != 3:
    print("Two arguments required!\n")
    quit()

try:
    old_file = open(sys.argv[1], "r")
    old_file.close()
    new_file = open(sys.argv[2], "r")
    new_file.close()
except:
    print("Need file permissions!\n")
    quit()

out_name = sys.argv[2] + ".csv"
out_file = open(out_name, "w")

old = pose_from_pdb(sys.argv[1])
new = pose_from_pdb(sys.argv[2])
num_res = old.total_residue()

out_file.write("Res #,Res Name,Original,Designed,Original # Adjacent,Designed # Adjacent\n")

old_bur_num = 0
old_bur_cons = 0
old_exp_num = 0
old_exp_cons = 0
for j in xrange(1, num_res + 1):
    old_res_CA = old.residue(j).xyz("CA")
    new_res_CA = new.residue(j).xyz("CA")
    old_num_adj = 0
    new_num_adj = 0
    for k in xrange(1, num_res + 1):
        if j != k and old.residue(k).name()[:3] != "GLY":
            old_res_CB = old.residue(k).xyz("CB")
            old_vector = old_res_CA - old_res_CB
            old_dist = old_vector.norm
            if old_dist <= 8:
                old_num_adj += 1
    for k in xrange(1, num_res + 1):
        if j != k and new.residue(k).name()[:3] != "GLY":
            new_res_CB = new.residue(k).xyz("CB")
            new_vector = new_res_CA - new_res_CB
            new_dist = new_vector.norm
            if new_dist <= 8:
                new_num_adj += 1
    if old_num_adj > 7 and new_num_adj > 7:
        out_file.write(str(j) + "," + old.residue(j).name()[:3] + ",B,B," + str(old_num_adj) + "," + str(new_num_adj) + "\n")
        old_bur_num += 1
        old_bur_cons += 1
    elif old_num_adj > 7 and new_num_adj <= 7:
        out_file.write(str(j) + "," + old.residue(j).name()[:3] + ",B,E," + str(old_num_adj) + "," + str(new_num_adj) + "\n")
        old_bur_num += 1
    elif old_num_adj <= 7 and new_num_adj > 7:
        out_file.write(str(j) + "," + old.residue(j).name()[:3] + ",E,B," + str(old_num_adj) + "," + str(new_num_adj) + "\n")
        old_exp_num += 1
    else:
        out_file.write(str(j) + "," + old.residue(j).name()[:3] + ",E,E," + str(old_num_adj) + "," + str(new_num_adj) + "\n")
        old_exp_num += 1
        old_exp_cons += 1

out_file.write(str(old_bur_cons) + "/" + str(old_bur_num) + " buried residues remained buried\n")
out_file.write(str(old_exp_cons) + "/" + str(old_exp_num) + " exposed residues remained exposed\n")
out_file.write("RMSD for alpha-C's is " + str(CA_rmsd(old, new)) + "\n")
out_file.write("RMSD over all atoms is " + str(all_atom_rmsd(old, new)) + "\n")
out_file.close()
