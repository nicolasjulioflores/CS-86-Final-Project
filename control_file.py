import time as time
from struc_predict_repack_singular import predict_repack

# start pyrosetta stuff
from rosetta import *
init()


if __name__ == "__main__":
  SEQ = "SWTWEGNKWTWK" # string sequence
  P_FILE = "1le0_m1.pdb" # pdb file
  OUTFILE = "rmsd_output_repack.txt"
  REPACK_FREQ = [10, 100, 1000, 10000]
  STANDARD = False

  # Write the header for the output file
  with open(OUTFILE, "w") as fout:
    fout.write("rmsd and score between current and native structure\nrmsd\t\tscore\n")

  # Write the header for the time file
  time_file = open("time_file_{0}.txt".format(P_FILE[:-4]), 'w')
  time_file.write(P_FILE + " standard: " + str(STANDARD) + "\n")

  for freq in REPACK_FREQ:
    start = time.time()
    predict_repack(SEQ, P_FILE, OUTFILE, freq, standard=STANDARD)
    end = time.time()

    real_time = end - start
    time_file.write(str(freq) + " " + str(real_time) + "\n") 
  
  time_file.close()