import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math

# Read in data
RMSDs = list()
with open("RMSDs.txt", "r") as f:
    for line in f:
        RMSDs.append(float(line.strip()))

scores = list()
with open("scores.txt", "r") as f:
    for line in f:
        scores.append(float(line.strip()))

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
plt.title('score3 vs RMSD for kT = 1')

plt.xlabel('RMSD')
plt.ylabel('score3')

plt.show()
