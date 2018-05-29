import time as time
from struc_predict_repack_singular import predict_repack

# start pyrosetta stuff
from rosetta import *
init()


if __name__ == "__main__":
  NUM_RUNS = 2
  SEQ = ["SWTWEGNKWTWK", "SWTWEGNKWTWK"] # string sequence
  P_FILE = ["1le0_m1.pdb", "1le0_m1.pdb"] # pdb file
  STANDARD = [True, False]
  OUTBASE = "rmsd_output_repack"
  REPACK_FREQ = [10, 100, 1000, 10000]
  
  for i in xrange(NUM_RUNS):
    
    if STANDARD[i]:
      out_s = "fullatom"
    else:
      out_s = "centroid"

    OUTFILE = OUTBASE + "_" + out_s + "_" + P_FILE[i][:-4] + ".txt"

    # Write the header for the output file
    with open(OUTFILE, "w") as fout:
      fout.write("rmsd and score between current and native structure\nrmsd\t\tscore\n")

    # Write the header for the time file
    time_file = open("time_file_{0}.txt".format(P_FILE[i][:-4]), 'w')
    time_file.write(P_FILE[i] + " standard: " + str(STANDARD[i]) + "\n")

    for freq in REPACK_FREQ:
      start = time.time()
      predict_repack(SEQ[i], P_FILE[i], OUTFILE, freq, standard=STANDARD[i])
      end = time.time()

      real_time = end - start
      time_file.write(str(freq) + " " + str(real_time) + "\n") 
    
    time_file.close()