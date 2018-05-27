import random
import math

# start pyrosetta stuff
from rosetta import *
init()


def predict_repack(in_sequence, in_protein, outfile, repack_freq, standard = False, length=100000, kT=1):
  
  # Full atom or centroid?
  if standard == False:
    method = 'centroid'
    scoring = 'score3'
  else:
    method = 'fa_standard'
    scoring = 'standard'

  # initial pose, based on sequence only and using given representation
  sim = pose_from_sequence(in_sequence, method)
  sim.pdb_info().name(inprotein[:-4] + "_sim")

  # actual 1le0 structure, taken from pymol
  correct = pose_from_pdb(inprotein)

  # create scorer using score3 and set an initial score
  scorer = create_score_function(scoring)
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
  bestScore = 10000000

  out_file = open(outfile, 'a')
  out_file.write("repacking frequency for this run is " + str(repack_freq[j]) + "\n")

  total_count = 0
  count = 0
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
    if (total_count % repack_freq == 0):
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
  dump_pdb(best, 'output_' + str(repack_freq) +'.pdb')

  out_file.close()