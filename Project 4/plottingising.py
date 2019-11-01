import numpy as np
import matplotlib.pyplot as plt

filename = "outfile.txt"
fil = open(filename, 'r') #reading file
lines = fil.readlines() #splitting in lines
length  = len(lines)


# Initialize lists
E = []
M = []
heatcapacity = []
T = []
X = []

for line in (lines):
    T.append(float(line.split()[0]))
    E.append(float(line.split()[1]))
    heatcapacity.append(float(line.split()[2]))
    X.append(float(line.split()[3]))
    M.append(float(line.split()[4]))
fil.close()

E = np.array(E)
T = np.array(T)
M = np.array(M)
C = np.array(heatcapacity)

plt.plot(T, E)
plt.title('Mean Energy')
plt.xlabel("Temperature")
plt.ylabel("Mean energy/n^2")
plt.grid()
plt.show()

plt.plot(T, M)
plt.title('Abs Mean Magnetization')
plt.xlabel("Temperature")
plt.ylabel("Mean absolute Magnetization/n^2")
plt.grid()
plt.show()

plt.plot(T, C)
plt.title('heatcapacity')
plt.xlabel("Temperature")
plt.ylabel("sigma_E ^2/ T^2")
plt.grid()
plt.show()