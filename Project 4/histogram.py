import numpy as np
import matplotlib.pyplot as plt

#compile >> g++ maind.cpp lib.cpp << and run
#this program. It will ask for the filename and temp to plot histogram.

filename = str(input("Filename : "))
fil = open(filename, 'r') #reading file
lines = fil.readlines() #splitting in lines
length  = len(lines)
T = float(input("T ="))

E = []

for line in lines:
    E.append(float(line.split()[0]))
fil.close()

E = np.array(E)
#each MC cycle represent a "time step"

#assuming that the energy is stabilized after 4000 MCc
#based on the plots from previou
E_stabilized = E[4000:-1]

#plotting
plt.figure(figsize = (8,6))

plt.title("Histogram showing energy probability P(E) per particle   \n \
20 x 20 lattice, T = %.1f   " % T,fontsize = 13)
plt.xlabel("Energy per particle E / N",fontsize = 13)
plt.ylabel("Probability P(E) / N", fontsize = 13)
plt.xticks(fontsize = 13)
plt.yticks(fontsize = 13)

bins = 30
weights = np.ones_like(E_stabilized)/float(len(E_stabilized))
plt.hist(E_stabilized, bins = 30, color = 'blue', edgecolor = 'black'\
,normed = 0, weights=weights)
#normalized histogram of energies

plt.show()

"""
sample run
Terminal > histogram.py
(output is a plot)
"""
