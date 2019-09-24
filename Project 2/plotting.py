import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
eigvects = np.zeros((N,N))
filename = "eigenvectors_%i.txt"%(N)
fil=open(filename, 'r') #reading file
line1 = fil.readline()
line2 = fil.readline()

lines=fil.readlines() #splitting in lines
N_new = len(lines)
if not N == N_new:
    print("Something went wrong, the eigenvectors are not of the correct size!")
    sys.exit(1)

for num,i in enumerate(lines): #putting everything in the right list
    text = i.split()
    eigvects[num,:] = [float(txt) for txt in text]
fil.close()

rho = np.linspace(0,1,N)
plt.plot(rho,eigvects[:,31])
plt.show()