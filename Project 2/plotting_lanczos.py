import numpy as np
import matplotlib.pyplot as plt
import sys

# Plots eigenvalue estimates and time used for Lanczos' algo
# Requires N as input, where the matrix one used Lanczos on is a NxN matrix

N = int(sys.argv[1])

fil = open("Lanczos_N_and_K.txt", "r")
lines = fil.readlines()
Ns = np.zeros(len(lines),dtype = "int")
ks = np.zeros_like(Ns)

#Read values for N and k
for i, line in enumerate(lines):
    text = line.split("=")
    _N = int(text[1].split()[0])
    k = int(text[2].split()[0])
    Ns[i] = _N; ks[i] = k

#Array that tells which lines are relevant for this plot
bool_array = Ns == N
fil.close()

plt.figure()
plt.title("Eigenvalue estimates, N = %i"%(N), fontsize = 16)
plt.xlabel("k", fontsize = 14)
plt.ylabel(r"$\lambda$", fontsize = 14)

eigvals_analytical = [3,7,11,15,18,22]

# Find and plot eigenvalues lower than 22 if it is from a NxN matrix
for i,k in enumerate(ks):
    if bool_array[i]:
        filename = "eigenvalues_%i_k%i.txt"%(N,k)
        fil = open(filename, 'r') #reading file
        fil.readline()
        fil.readline()
        eigvals = []
        line = fil.readline()
        eigval = float(line)
        while eigval<=22:
            eigvals.append(eigval)
            line = fil.readline()
            eigval = float(line)
        plt.semilogx(np.ones(len(eigvals))*k, eigvals, '+',markersize=10,linewidth=4)

#Plot analytical eigenvalues
for eigval in eigvals_analytical:
    plt.semilogx([min(ks),max(ks)],[eigval,eigval],'k--')

plt.savefig("lanczos_eigvals_%i.png"%(N))

plt.figure()
plt.title("Time used, Lanczos' and QR, N = %i"%(N), fontsize = 16)
plt.xlabel("k", fontsize = 14)
plt.ylabel("t [s]", fontsize = 14)

#Manually copy-pasted results from terminal
k = [50,60,70,80,90,100,150,200,300,400,500,600,700,800,900,1000]
t = [1.96875,2.32812,2.71875,3.26562,3.60938,3.89062,5.84375,9.34375,18.0469,28.7031,48.6562,75.9219,56.7031,75.1875,99.5312,125.969]
plt.loglog(k,t,'x',markersize=8)
plt.savefig("Lanczos_time_N%i.png"%(N))

plt.show()