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

    for L in np.arange(40,101,20):
        filename = f"outfile_MPI_L{L}_n1e6.txt"
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

    plt.figure('E')
    plt.legend()

    plt.figure('Mag')
    plt.legend()

    plt.figure('hc')
    plt.legend()

    plt.show()

load_plot_mpi()

"""
Plotting of the energy and magnetization over time
"""
filename = "temp.txt"
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
x = np.arange(1, 10000, 10)
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
plt.show()