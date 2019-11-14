import numpy as np
import matplotlib.pyplot as plt

filename = "histdata.txt"
fil = open(filename, 'r') #reading file
lines = fil.readlines() #splitting in lines
length  = len(lines)
T = float(input("T ="))


MCc = []
E = []
for line in lines:
    MCc.append(float(line.split()[0]))
    E.append(float(line.split()[1]))
fil.close()

#each MC cycle represent a "time step"
E = np.array(E)
MCc = np.array(MCc)

plt.figure(figsize = (8,6))
plt.title("Energy plotted against MC cycles  \n \
20 x 20 lattice, T = %.1f   " % T,fontsize = 12)
plt.xlabel("MC cycles",fontsize = 12)
plt.ylabel("Energy per particle E / N", fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.plot(MCc[0::1500],E[0::1500])


#assuming that the energy is stabilized after 3000 MCc
E_stabilized = E[3000:-1]
plt.figure(figsize = (8,6))
plt.title("Histogram showing Probability density P(E) per particle   \n \
20 x 20 lattice, T = %.1f   " % T,fontsize = 12)
plt.xlabel("Energy per particle E / N",fontsize = 12)
plt.ylabel("Probability density P(E) / N", fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)

weights = np.ones_like(E_stabilized)/float(len(E_stabilized))
plt.hist(E_stabilized, bins = 30, color = 'blue', edgecolor = 'black',normed = 0, weights=weights)
#normalized histogram
plt.show()
