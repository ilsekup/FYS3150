import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
eigvects = np.zeros((N,N))
filename = "eigenvectors_%i.txt"%(N)
fil = open(filename, 'r') #reading file
fil.readline()
fil.readline()

lines = fil.readlines() #splitting in lines
N_new = len(lines)
if not N == N_new:
    print("Something went wrong, the eigenvectors are not of the correct size!")
    sys.exit(1)

for num,i in enumerate(lines): #Putting the eigenvectors into arrays
    text = i.split()
    eigvects[num,:] = [float(txt) for txt in text]
fil.close()

eigenvalues = np.zeros(N)
filename = "eigenvalues_%i.txt"%(N)
fil=open(filename,'r')
fil.readline()
fil.readline()
lines = fil.readlines()
for num,i in enumerate(lines): #putting the eigenvalues into lists
    eigenvalues[num] = float(i)
fil.close()

rho = np.linspace(0,10,N)

#Plots the eigenvectors for these eigenvalues (indexed by number)
n_plot = [0,1,2]

plt.figure(figsize = (7,5))
plt.xlabel(r"$\rho$", fontsize = 14)
plt.ylabel(r"$u(\rho)$", fontsize = 14)

for n in n_plot:
    plt.plot(rho,eigvects[:,n],label = r"$\lambda_%i \approx %.2f$"%(n,eigenvalues[n]))
plt.legend(fontsize = 12)
plt.savefig("solutions.png")
plt.show()