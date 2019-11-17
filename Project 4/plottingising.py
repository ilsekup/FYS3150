import numpy as np
import matplotlib.pyplot as plt


def load_plot_mpi():
    """
    loads data and plots data from mpi_ising program from run_mpi.sh
    """
    plt.figure('E')
    plt.title('Mean Energy')
    plt.xlabel("Temperature")
    plt.ylabel("Mean energy/n^2")
    plt.grid()

    plt.figure('Mag')
    plt.title('Abs Mean Magnetization')
    plt.xlabel("Temperature")
    plt.ylabel("Mean absolute Magnetization/n^2")
    plt.grid()


    plt.figure('hc')
    plt.title('heatcapacity')
    plt.xlabel("Temperature")
    plt.ylabel(r"$\sigma_E^2/ T^2$")
    plt.grid()

    plt.figure('X')
    plt.title('Suceptibility')
    plt.xlabel("Temperature")
    plt.ylabel(r"$\sigma_M^2/T^2$")
    plt.grid()

    for L in np.arange(40,101,20):
        filename = f"outfile_MPI_L{L}.txt"
        fil = open(filename, 'r') #reading file
        lines = fil.readlines()
        l = len(lines)
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
        E = np.array(E)
        T = np.array(T)
        M = np.array(M)
        C = np.array(heatcapacity)
        X = np.array(X)
        fil.close()
        plt.figure('E')
        plt.plot(T, E,label=f"L = {L}")

        plt.figure('Mag')
        plt.plot(T, M,label=f"L = {L}")
    
        plt.figure('hc')
        plt.plot(T, C,label=f"L = {L}")

        plt.figure('X')
        plt.plot(T, X,label=f"L = {L}")

    plt.figure('E')
    plt.legend()
    plt.savefig("Energy.pdf")

    plt.figure('Mag')
    plt.legend()
    plt.savefig("Magnetism.pdf")

    plt.figure('hc')
    plt.legend()
    plt.savefig("Heatcapacity.pdf")

    plt.figure('X')
    plt.legend()
    plt.savefig("Susceptibility.pdf")
        
    plt.show()

load_plot_mpi()
"""
filename = "outfilempi.txt"
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
X = np.array(X)

plt.figure()
plt.plot(T, E)
plt.title('Mean Energy')
plt.xlabel("Temperature")
plt.ylabel("Mean energy/n^2")
plt.grid()

plt.figure()
plt.plot(T, M)
plt.title('Abs Mean Magnetization')
plt.xlabel("Temperature")
plt.ylabel("Mean absolute Magnetization/n^2")
plt.grid()

plt.figure()
plt.plot(T, C)
plt.title('heatcapacity')
plt.xlabel("Temperature")
plt.ylabel("sigma_E ^2/ T^2")
plt.grid()
plt.show()
"""
# mc =  [500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 30000, 40000,50000, 100000]
"""
mc = np.arange(0, 100000, step=10000)
Etemp1 = np.zeros(len(mc)); Mtemp1 =np.zeros(len(mc)); Etemp24 = np.zeros(len(mc)); Mtemp24 =np.zeros(len(mc))
Etemp1ran = np.zeros(len(mc)); Mtemp1ran =np.zeros(len(mc))
"""
# for i, s in enumerate(mc):
#     filename = "outfile%g.txt" %(s)
#     fil = open(filename, 'r') #reading file
#     lines = fil.readlines() #splitting in lines
#     length  = len(lines)
#     E = []
#     M = []
#     for line in (lines):
#         E.append(float(line.split()[1]))
#         M.append(float(line.split()[-1]))
#     fil.close()
#     E = np.array(E)
#     M = np.array(M)
#     Etemp1[i]  = E[0];
#     Etemp24[i] = E[14]
#     Mtemp1[i] = M[0];
#     Mtemp24[i] = M[14]
"""
for i, s in enumerate(mc):
    filename = "outfileran%g.txt" %(s)
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
    Etemp1ran[i]  = E[0]
    Mtemp1ran[i] = M[0]
mc = np.array(mc)
fig, axs = plt.subplots(2, sharex=True)
fig.suptitle('Temperature T=1.0 ')
axs[0].plot(x, E, label='ordered')
axs[0].set_title('Energy')
axs[0].legend()
axs[0].set_ylabel('Energy/J')
axs[1].plot(x, M, label='ordered')
axs[1].set_title('Magnetization')
axs[1].legend()
axs[1].set_ylabel('|M|')
plt.xlabel('Monte Carlo cycles/ time')
# plt.xticks(mc)
# plt.xscale('log')
#plt.ticklabel_format(useOffset=False, style='plain')
plt.show()
# plt.plot(Etemp1, np.linspace(0, 19, 19))
# plt.show()
"""
