import numpy as np
import matplotlib.pyplot as plt
import sys

error = np.zeros(9)
integrationnumber = [20, 30, 50, 70, 100, 150, 200, 250, 300] #integrationnumber I will look at
for i, s in enumerate(integrationnumber):
    filename = "eigenvalues_%i.txt"%(s)
    fil=open(filename, 'r') #reading file
    line1 = fil.readline()
    error[i] = float(line1.split(" ")[-1])
    lines = fil.readlines()
error_nopot = np.zeros(9)
for i, s in enumerate(integrationnumber):
    filename = "eigenvalues_nopot%i.txt"%(s)
    fil=open(filename, 'r') #reading file
    line1 = fil.readline()
    error_nopot[i] = float(line1.split(" ")[-1])
    lines = fil.readlines()

error_rho100 = np.zeros(9)
for i, s in enumerate(integrationnumber):
    filename = "eigenvalues_rho100%i.txt"%(s)
    fil=open(filename, 'r') #reading file
    line1 = fil.readline()
    error_rho100[i] = float(line1.split(" ")[-1])
    lines = fil.readlines()

error_rho4 = np.zeros(9)
for i, s in enumerate(integrationnumber):
    filename = "eigenvalues_rho4%i.txt"%(s)
    fil=open(filename, 'r') #reading file
    line1 = fil.readline()
    error_rho4[i] = float(line1.split(" ")[-1])
    lines = fil.readlines()

integrationnumber = np.array(integrationnumber)
plt.semilogy(integrationnumber, error, label='with potential inf = 10', color='b')
plt.semilogy(integrationnumber, error_rho4, label='with potential inf = 4', color='r')
plt.semilogy(integrationnumber, error_rho100, label='with potential inf = 100', color='k')
plt.title('Error for different infinity approximations')
plt.xlabel('Integrationspoints N')
plt.ylabel('Log of error')
plt.legend()
plt.grid()
plt.show()

#plt.plot(integrationnumber, error_nopot, label='without potential', color='c')

needed = np.array([437, 1163, 3707, 7433, 15578, 35816, 63930,  100393, 145504])
needed_nopot = np.array([655, 1535, 4349, 8621, 17660, 39890, 70833, 110976, 159933])
plt.plot(integrationnumber, needed, label='With potential inf = 10', color='g')
plt.plot(integrationnumber, needed_nopot, label='Without potential', color='r')
plt.title('Number of Iterations Needed')
plt.xlabel('Integrationspoints N')
plt.legend()
plt.grid()
plt.show()

#semilogy