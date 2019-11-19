import numpy as np
import matplotlib.pyplot as plt

filename = "mcsweeps.txt"
fil = open(filename, 'r') #reading file
lines = fil.readlines() #splitting in lines
length  = len(lines)
E = []
M = []
sweeps = []
counts = []
for line in (lines):
    sweeps.append(float(line.split()[0]))
    E.append(float(line.split()[1]))
    M.append(float(line.split()[2]))
    counts.append(float(line.split()[3]))
fil.close()

E = np.array(E)
M = np.array(M)
sweeps = np.array(sweeps)
counts = np.array(counts)

T = float(input("T = "))

#Plotting of the energy and magnetization over time
fig, axs = plt.subplots(2, sharex=True)
fig.suptitle('Temperature T= %.1f' % T)
axs[0].plot(sweeps, E, label='ordered')
axs[0].set_title('Energy')
axs[0].legend()
axs[0].set_ylabel('Energy/J')
axs[1].plot(sweeps, M, label='ordered')
axs[1].set_title('Magnetization')
axs[1].legend()
axs[1].set_ylabel('|M|')
plt.xlabel('Monte Carlo cycles/ time')

#Plotting accepct configs. as a function a MC sweeps
plt.figure()
plt.title("Number of accepted moves vs MC sweeps \n 20 x 20 lattice at T = %.1f" % T, \
fontsize = 12)
plt.plot(sweeps,counts)
plt.xlabel("Number of MC sweeps", fontsize = 12)
plt.ylabel("Number of accepted configurations", fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.show()
