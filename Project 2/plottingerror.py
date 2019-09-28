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
integrationnumber = [20, 30, 50, 70, 100, 150, 200, 250, 300] #integrationnumber I will look at
for i, s in enumerate(integrationnumber):
    filename = "eigenvalues_nopot%i.txt"%(s)
    fil=open(filename, 'r') #reading file
    line1 = fil.readline()
    error_nopot[i] = float(line1.split(" ")[-1])
    lines = fil.readlines()


integrationnumber = np.array(integrationnumber)
plt.plot(integrationnumber, error, label='with potential', color='b')
plt.plot(integrationnumber, error_nopot, label='without potential', color='c')
plt.title('Error for the first 20 values')
plt.xlabel('Integrationspoints N')
plt.ylabel('Relative Error between analytical and calculated')
plt.legend()
plt.show()

needed = np.array([437, 1163, 3707, 7433, 15578, 35816, 63930,  100393, 145504])
needed_nopot = np.array([655, 1535, 4349, 8621, 17660, 39890, 70833, 110976, 159933])
plt.plot(integrationnumber, needed, label='With potential', color='g')
plt.plot(integrationnumber, needed_nopot, label='Without potential', color='r')
plt.title('Number of Iterations Needed')
plt.xlabel('Integrationspoints N')
plt.legend()
plt.show()
