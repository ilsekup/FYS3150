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

# plt.plot(T, E)
# plt.title('Mean Energy')
# plt.xlabel("Temperature")
# plt.ylabel("Mean energy/n^2")
# plt.grid()
# plt.show()
#
# plt.plot(T, M)
# plt.title('Abs Mean Magnetization')
# plt.xlabel("Temperature")
# plt.ylabel("Mean absolute Magnetization/n^2")
# plt.grid()
# plt.show()
#
# plt.plot(T, C)
# plt.title('heatcapacity')
# plt.xlabel("Temperature")
# plt.ylabel("sigma_E ^2/ T^2")
# plt.grid()
# plt.show()

mc = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
Etemp1 = np.zeros(len(mc)); Mtemp1 =np.zeros(len(mc))
for i, s in enumerate(mc):
    filename = "outfile%g.txt" %(s)
    fil = open(filename, 'r') #reading file
    lines = fil.readlines() #splitting in lines
    length  = len(lines)
    E = []
    M = []
    for line in (lines):
        E.append(float(line.split()[1]))
        M.append(float(line.split()[-1]))
    fil.close()
    E = np.array(E)
    M = np.array(M)
    Etemp1[i]  = E[0]
    Mtemp1[i] = M[0]
mc = np.array(mc)
fig, axs = plt.subplots(2, sharex=True)
fig.suptitle('Temperature 1 ')
axs[0].plot(mc, Etemp1)
axs[0].set_title('Energy')
axs[1].plot(mc, Mtemp1)
axs[1].set_title('Magnetization')
plt.xticks(mc)
plt.xscale('log')

#plt.ticklabel_format(useOffset=False, style='plain')
plt.show()