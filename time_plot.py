import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from scipy import stats

import numpy as np

TIME_DIR = "time/"
OUTPUT_DIR = "output/"
NORM_FREQ = 10000
NORMALIZE = False
PROTEINS = {'1le0':'blue', '1ubq':'cyan', '1x8y':'darkred', '5w0h':'red', '6elm':'olive', '4xdx':'green', '4ytv':'purple', '4yud':'yellow'}


def read_time_file(infile, normalize):
  data = list()

  with open(infile, 'r') as f:
    lines = f.readlines()
  
  
  for i in xrange(1, len(lines)):
    pair = lines[i].split()

    repack_freq = float(pair[0].strip())
    time = float(pair[1].strip())

    data.append(time)

    if repack_freq == NORM_FREQ:
      norm = time

  if normalize:
    # Normalize the data
    data = [float(x)/float(norm) for x in data]

  return data

def read_output_file(outfile, normalize):
  '''
  Purpose: read in the output files
  '''
  
  def read_to_dict(outfile, dictionary):
    with open(outfile, 'r') as f:

      first = True
      for line in f:
        if not line.startswith("rmsd"):
          if line.startswith("repacking"):
            
            if not first:
              dictionary['RMSDs'].append(RMSDs)
              dictionary['scores'].append(scores)
              dictionary['best_scores'].append(min(scores))
              dictionary['best_RMSD'].append(min(RMSDs))
            RMSDs = list()
            scores = list()
            first = False

          else:
            pair = line.split()
            
            try: 
              float(pair[0].strip())
              float(pair[1].strip())
            except ValueError:
              continue

            RMSDs.append(float(pair[0].strip()))
            scores.append(float(pair[1].strip()))
      
    dictionary['RMSDs'].append(RMSDs)
    dictionary['scores'].append(scores)
    dictionary['best_scores'].append(min(scores))
    dictionary['best_RMSD'].append(min(RMSDs))

    if normalize:
      norm_1 = min(dictionary['best_scores'])
      dictionary['best_scores'] = [x + norm_1 + 1 for x in dictionary['best_scores']]

      norm_2 = min(dictionary['best_RMSD'])
      dictionary['best_RMSD'] = [float(x)/float(norm_2) for x in dictionary['best_RMSD']]
    
    return dictionary

  method = dict()
  method['RMSDs'] = list()
  method['scores'] = list()
  method['best_RMSD'] = list()
  method['best_scores'] = list()

  return read_to_dict(outfile, method)


def _generate_x_y(time_data_dictionary):
  
  X = list()
  Y = list()

  for key in time_data_dictionary.keys():
    for value in time_data_dictionary[key]:
      X.append(key)
      Y.append(value)
  
  return X, Y

def plot_time_data(time_dir, output_dir, normalize=False):
  
  # Needed for plotting
  frequencies = [10, 100, 1000, 10000]

  centroid = dict()
  fullatom = dict()

  for item in [centroid, fullatom]:
    item['RMSD'] = dict()
    item['Score'] = dict()
    item['Time'] = dict()

  # Read in time data
  for file in os.listdir(time_dir):
    datum = read_time_file(os.path.join(time_dir, file), normalize)

    for name in PROTEINS.keys():
      if name in file and datum != None:
        
        if "fullatom" in file:
          fullatom['Time'][name] = datum
        elif "centroid" in file:
          centroid['Time'][name] = datum

  # Read in output data
  for file in os.listdir(output_dir):
    datum = read_output_file(os.path.join(output_dir, file), normalize)
    
    for name in PROTEINS.keys():
      if name in file and datum != None:

        if "fullatom" in file:
          fullatom['RMSD'][name] = datum['best_RMSD']
          fullatom['Score'][name] = datum['best_scores']
        elif "centroid" in file:
          centroid['RMSD'][name] = datum['best_RMSD']
          centroid['Score'][name] = datum['best_scores']

  # Plot the data for each method
  avg_table = dict()
  for method in [fullatom, centroid]:
    if method == fullatom:
      method_name = 'Fullatom'
    else:
      method_name = 'Centroid'

    avg_table[method_name] = dict()
    for key in method.keys():
      
      avg_table[method_name][key] = dict()
      patches = list()
      x_total = list()
      y_total = list()
      for name in PROTEINS.keys():
        if name in method[key].keys(): 
          x = frequencies
          y = method[key][name]

          plt.scatter(x, y)

          x_total.append(x)
          y_total.append(y)

          patches.append(mpatches.Patch(color=PROTEINS[name], label=name))
      
      x_total = np.array(x_total).flatten()
      y_total = np.array(y_total).flatten()
      
      for freq in frequencies:
        y_freq = [y for x,y in zip(x_total, y_total) if x == freq]

        avg_table[method_name][key][freq] = [np.mean(y_freq), np.std(y_freq)]

      line = np.polyfit(x_total, y_total, 2)
      line_fn = np.poly1d(line)

      plt.plot([], [],'o',frequencies, line_fn(frequencies),'k')

      plt.title('{0} scoring: {1} vs Repacking Frequency'.format(method_name, key))
      plt.xscale('log')
      plt.xlabel('Repacking Frequency')
      plt.ylabel(key)
      plt.legend(handles=patches)

      plt.show()
      plt.clf()

    print("Table for: {}".format(method_name))
    print("")
    print("\\begin{tabular}{| p{2cm}| c | c | c | c |}")
    print("\\hline")
    print("Repacking Frequency $& 10^1 & 10^2 & 10^3 & 10^4$\\\\")
    print("\\hline")

    # Rounding function used to round the raw floats
    for y_type in avg_table[method_name].keys():
      total_string = "{0} & ".format(y_type)
      for freq in frequencies:
        data = avg_table[method_name][y_type][freq]
        data_string = "{0} $\\pm$ {1} & ".format(round(data[0],3), round(data[1],3))
        total_string += data_string

      print(total_string[:-2] + " \\\\")

    print("\\hline")
    print("\\end{tabular}")
    print("")

  

    # # Plot times
    # patches = list()
    # for name in PROTEINS.keys():
    #   if name in method['time'].keys(): 
    #     x_1 = frequencies
    #     y_1 = method['time'][name] 
    #     plt.scatter(x_1, y_1)
        
    #     patches.append(mpatches.Patch(color=PROTEINS[name], label=name))
    

    # plt.title('{0} scoring: Time vs Repacking Frequency'.format(method_name))
    # plt.xscale('log')
    # plt.xlabel('Repacking Frequency')
    # plt.ylabel('Time')
    # plt.legend(handles=patches)
    # plt.show()
    # plt.clf()

    # # Plot scores
    # patches = list()
    # for name in PROTEINS.keys():
    #   if name in method['scores'].keys():
    #     x_2 = frequencies
    #     y_2 = method['scores'][name]
    #     plt.scatter(x_2, y_2, color=PROTEINS[name])
        
    #     patches.append(mpatches.Patch(color=PROTEINS[name], label=name))

    # plt.title("{0} scoring: Score vs Repacking Frequency".format(method_name))
    # plt.xscale('log')
    # plt.xlabel('Repacking Frequency')
    # plt.ylabel('Score')
    # plt.legend(handles=patches)
    # plt.show()
    # plt.clf()

    # # Plot RMSD
    # patches = list()
    # for name in PROTEINS.keys():
    #   if name in method['RMSD'].keys():
    #     x_3 = frequencies
    #     y_3 = method['RMSD'][name]

    #     plt.scatter(x_3, y_3, color=PROTEINS[name])
    #     patches.append(mpatches.Patch(color=PROTEINS[name], label=name))
   
    # plt.title("{0} scoring: RMSD vs Repacking Frequency".format(method_name))
    # plt.xscale('log')
    # plt.xlabel('Repacking Frequency')
    # plt.ylabel('RMSD')
    # plt.legend(handles=patches)
    # plt.show()
    # plt.clf()


if __name__ == "__main__":
  plot_time_data(TIME_DIR, OUTPUT_DIR, NORMALIZE)

  
