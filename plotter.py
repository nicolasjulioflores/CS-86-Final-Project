import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math

OUTFILE = 'rmsd_output_repack.txt'

def read_output_file(outfile):
    '''
    Purpose: read in the output files
    '''

    RMSDs = list()
    scores = list()
    all_RMSD = list()
    all_scores = list()
    with open(outfile, 'r') as f:
        for line in f:
            if not line.startswith("rmsd"):
                if line.startswith("repacking"):
                    all_RMSD.append(RMSDs)
                    all_scores.append(scores)
                    RMSDs = list()
                    scores = list()
                else:
                    pair = line.split()
                    RMSDs.append(float(pair[0].strip()))
                    scores.append(float(pair[1].strip()))
    
    all_RMSD.append(RMSDs)
    all_scores.append(scores)

    return all_RMSD, all_scores
            
def plot_RMSD_vs_score(RMSDs, scores, title_string=None):
    '''
    Purpose: plot the RMSD vs score for a given run
    '''

    # Colors to be used in the plot from left to right
    colors = ['darkred', 'red', 'orangered', 'orange', 'gold', 'yellow']

    # Create bins for plotting
    n = len(scores)
    bins = [int(math.floor(i)) for i in np.linspace(0, n, len(colors) + 1)]

    # Values are plotted by color according to their step #
    for i in xrange(len(bins) - 1):
        plt.scatter(RMSDs[bins[i]:bins[i+1]], scores[bins[i]:bins[i+1]],
                    c=colors[i]) 

    # Last value is plotted as a blue dot
    plt.plot(RMSDs[-1], scores[-1], 'bo')

    # Legend created using: https://matplotlib.org/users/legend_guide.html
    patches = [mpatches.Patch(color=colors[i], label=str(bins[i]) + " to " 
                + str(bins[i+1])) for i in xrange(len(bins) - 1)]
    last_one = [mpatches.Patch(color='blue', label='last recorded value')]

    plt.legend(handles=patches+last_one)


    # Axes + Figure titles
    if title_string == None:
        plt.title('score vs RMSD')
    else 
        plt.title(title_string)

    plt.xlabel('RMSD')
    plt.ylabel('score')

    plt.show()

if __name__ == "__main__":
    
    # Read in the data
    all_RMSD, all_scores = read_output_file(OUTFILE)

    # The first entry in the list is always empty
    for i in xrange(1, len(all_RMSD)):
        plot_RMSD_vs_score(all_RMSD[i], all_scores[i])
        