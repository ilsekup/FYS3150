import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
eigvects = np.zeros((N,N))
filename = "eigenvectors_%i.txt"%(N)
fil = open(filename, 'r') #reading file

for i in range(2):
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

for i in range(3):
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
plt.ylabel(r"$P(\rho)$", fontsize = 14)
plt.title("The first %i solutions with N = %i points"%(len(n_plot),N),fontsize = 16)

for n in n_plot:
    plt.plot(rho,eigvects[:,n]**2,label = r"$\lambda_%i \approx %.2f$"%(n+1,eigenvalues[n]))

plt.legend(fontsize = 12)
plt.axis([0,4.5,0,0.022])
plt.savefig("solutions_N%i.png"%(N))
plt.show()