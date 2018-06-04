import os
import matplotlib.pyplot as plt

TIME_DIR = "time/"
OUTPUT_DIR = "output/"
NORM_FREQ = 10000
NORMALIZE = False

def read_time_file(infile, normalize):
  data = dict()

  with open(infile, 'r') as f:
    lines = f.readlines()
  
  
  for i in xrange(1, len(lines)):
    pair = lines[i].split()

    repack_freq = int(pair[0].strip())
    time = float(pair[1].strip())

    data[repack_freq] = time

    if repack_freq == NORM_FREQ:
      norm = time

  if normalize:
    # Normalize the data
    for key in data.keys():
      data[key] /= norm

  return data

def read_output_file(outfile, normalize):
  '''
  Purpose: read in the output files
  '''
  RMSDs = list()
  scores = list()
  all_RMSD = list()
  all_scores = list()
  best_scores = list()
  best_RMSD = list()
  with open(outfile, 'r') as f:
    for line in f:
      if not line.startswith("rmsd"):
        if line.startswith("repacking"):
          all_RMSD.append(RMSDs)
          all_scores.append(scores)
          if len(scores) != 0:
            best_scores.append(min(scores))
            best_RMSD.append(min(RMSDs))
          RMSDs = list()
          scores = list()
        else:
          pair = line.split()
          
          try: 
            float(pair[0].strip())
            float(pair[1].strip())
          except ValueError:
            continue

          RMSDs.append(float(pair[0].strip()))
          scores.append(float(pair[1].strip()))
    
  all_RMSD.append(RMSDs)
  all_scores.append(scores)
  best_scores.append(min(scores))
  best_RMSD.append(min(RMSDs))

  if normalize:
    norm_1 = min(best_scores)
    best_scores = [float(x)/float(norm_1) for x in best_scores]

    norm_2 = min(best_RMSD)
    best_RMSD = [float(x)/float(norm_2) for x in best_RMSD]

  return all_RMSD[1:], all_scores[1:], best_scores, best_RMSD

def _generate_x_y(time_data_dictionary):
  
  X = list()
  Y = list()

  for key in time_data_dictionary.keys():
    for value in time_data_dictionary[key]:
      X.append(key)
      Y.append(value)
  
  return X, Y

def plot_time_data(time_dir, output_dir, normalize=False):
  
  # Read in time data
  time_data = dict()
  for file in os.listdir(time_dir):
    datum = read_time_file(os.path.join(time_dir, file), normalize)

    for key, value in datum.iteritems():
      if key not in time_data.keys():
        time_data[key] = list()
      
      time_data[key].append(value)

  x_1, y_1 = _generate_x_y(time_data)

  # Read in output data
  frequencies = [10, 100, 1000, 10000]
  output_data = dict()
  RMSD_data = dict()
  for file in os.listdir(output_dir):
    a,b,output_datum, RMSD_datum = read_output_file(os.path.join(output_dir, file), normalize)
    
    for score, RMSD, freq in zip(output_datum, RMSD_datum, frequencies):
      if freq not in output_data.keys():
        output_data[freq] = list()
        RMSD_data[freq] = list()
      
      output_data[freq].append(score)
      RMSD_data[freq].append(RMSD)

  x_2, y_2 = _generate_x_y(output_data)

  x_3, y_3 = _generate_x_y(RMSD_data)

  plt.scatter(x_1, y_1)
  plt.xscale('log')
  plt.xlabel('Repacking Frequency')
  plt.ylabel('Time')
  plt.show()
  plt.clf()

  plt.scatter(x_2, y_2, color='red')
  plt.xscale('log')
  plt.xlabel('Repacking Frequency')
  plt.ylabel('Score')
  plt.show()
  plt.clf()

  plt.scatter(x_3, y_3, color='green')
  plt.xscale('log')
  plt.xlabel('Repacking Frequency')
  plt.ylabel('RMSD')
  plt.show()

if __name__ == "__main__":
  plot_time_data(TIME_DIR, OUTPUT_DIR, NORMALIZE)

  
