import random
import math

# start pyrosetta stuff
from rosetta import *
init()

# pymol = PyMOL_Mover()
# pymol.link.udp_ip = '10.31.112.1'   # IP address of my own machine

# initial pose, based on sequence only and using cetroid representation
sim = pose_from_sequence('SWTWEGNKWTWK','centroid')
sim.pdb_info().name('1le0_sim')
# pymol.apply(sim)

# actual 1le0 structure, taken from pymol
correct = pose_from_pdb('1le0_m1.pdb')

# create scorer using score3 and set an initial score
scorer = create_score_function('score3')
currentScore = scorer(sim)

# set up a temp pose 
temp = Pose()
temp.assign(sim)

# set up repacking task
task = standard_packer_task(temp)
task.restrict_to_repacking() # should we restrict or not?
mover = PackRotamersMover(scorer, task)

# set up a best pose
best = Pose()

# initial kT set to 1
kT = 1
# want 1000000 steps of MC
length = 1000000
# set a very high best score to initally compare against
bestScore = 1000000
# repacking frequency (repack every how many step)
repack_freq = [10, 100, 1000, 10000, 10000]

count = 0
total_count = 0

# create file to output rmsd and score values to for every 100th accepted structure
out_file = open('rmsd_output_repack.txt', 'w')
out_file.write("rmsd and score between current and native structure\nrmsd\t\tscore\n")

print "\nNOW STARTING MC SEARCH\n"

# go through all give repacking frequences

for j in range(0, len(repack_freq)):
  out_file.write("repacking frequency for this run is " + str(repack_freq[j]) + "\n")

  for i in range(length): 
    # find random residue number to change
    index = random.randint(1, 12) 

    #changing phi and psi angles in range [-10, 10]
    temp.set_phi(index,temp.phi(index) + random.randrange(-10,10))
    temp.set_psi(index,temp.psi(index) + random.randrange(-10,10))

    # score new conformation
    tempScore = scorer(temp)

    # check acceptence criteria
    # either temp score is better, or accept with defined probability
    if ((tempScore < currentScore) or (random.random() < math.exp(-(tempScore - currentScore)/kT))):
      # change kT = kT/1.05
      kT = kT/1.05

      # if temp score better than best score up to now, update best
      if (tempScore < bestScore):
        bestScore=tempScore
        best.assign(temp)
      
      # increment count of accepted structures
      count += 1

      # update current score and structure
      currentScore = tempScore
      sim.assign(temp)

      # pymol.apply(sim)

      #print out rmsd and score for every 100th accepted structure
      if (count%100 == 0):
        # get rmsd, write out rmsd and score
        dist = CA_rmsd(correct, sim)
        out_file.write(str(dist) + " ")
        out_file.write(str(currentScore) + "\n")

    # if not accepted, update temp
    else:
      # change kT = kT*1.01
      kT = kT * 1.01
      temp.assign(sim)

    total_count = total_count + 1
    if (total_count%repack_freq[j] == 0):
      mover.apply(temp)

  # inform about best score after MC run
  print "best score this run is " + str(scorer(best)) + "\n"
  print "it has a RMSD of " + str(CA_rmsd(correct, best))

  # write out best score to file
  dist = CA_rmsd(correct, best) 
  out_file.write(str(dist) + " ")
  out_file.write(str(scorer(best)) + "\n")

  print "\nSaving best structure of this run to pdb file\n"

  # put best pose into a new pdb file
  dump_pdb(best, 'output_' + str(repack_freq[j]) +'.pdb')

out_file.close()