import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def plot_mpi():
    """
    loads data and plots data from mpi_ising program from run_mpi.sh
    Also calculates T_C, and compares it to the actual value
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

    i = 0
    max_T = []
    max_T_hc = []
    L_ = []
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
        T_cubic = np.linspace(np.min(T),np.max(T),10000)
        fil.close()
        plt.figure('E')
        plt.plot(T, E,'x',label=f"L = {L}")

        plt.figure('Mag')
        plt.plot(T, M,'x',label=f"L = {L}")

        plt.figure('hc')
        fhc = interp1d(T,C, kind = 'cubic')
        plt.plot(T, C,f"xC{i}", label=f"L = {L}")
        plt.plot(T_cubic,fhc(T_cubic),f'--C{i}')
        max_T_hc.append(T_cubic[np.argmax(fhc(T_cubic))])

        plt.figure('X')
        # Using cubic splines to get a better estimate of T_C
        fX = interp1d(T, X, kind='cubic')
        plt.plot(T, X,f"xC{i}", label=f"L = {L}")
        plt.plot(T_cubic,fX(T_cubic),f'--C{i}')
        max_T.append(T_cubic[np.argmax(fX(T_cubic))])
        L_.append(L)
        i+=1

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

    plt.figure()
    L_inv = 1/np.array(L_)
    #Fit 1/T to T_C
    p,covmat = np.polyfit(L_inv,max_T,1,full = False, cov = True)
    p_hc,covmat_hc = np.polyfit(L_inv,max_T_hc,1,full = False, cov = True)

    plt.plot(L_inv,max_T,'xb', label = r"Data from $\chi$")
    plt.plot(L_inv,max_T_hc,'xr', label = "Data from heatcapacity")

    points = np.array([0,np.max(L_inv)])
    line = np.polyval(p,points)
    line_hc = np.polyval(p_hc,points)

    plt.plot(points,line,'--b')
    plt.plot(points,line_hc,'--r')

    b = p[-1]
    b_hc = p_hc[-1]

    sigma_b = np.sqrt(covmat[-1,-1])
    sigma_b_hc = np.sqrt(covmat_hc[-1,-1])

    print(f"Line intercept from chi = {b} ± {sigma_b}")
    print(f"Line intercept from heatcapacity = {b_hc} ± {sigma_b_hc}")
    act_val = 2/np.log(1+np.sqrt(2))
    sig_dist = np.abs(b-act_val)/sigma_b
    sig_dist_hc = np.abs(b_hc-act_val)/sigma_b_hc
    plt.xlabel("1/L")
    plt.ylabel(r"$T_C$")
    print(f"{sig_dist:.1f} sigma from the true value (chi)")
    print(f"{sig_dist_hc:.1f} sigma from the true value (heatcapacity)")

    plt.legend()
    plt.grid()
    plt.savefig("TC.pdf")

    plt.show()

def plot_mpi_timing():
    Ls = [40, 100]
    nps = [1,2,3,4]
    fix, axs = plt.subplots(2,2)

    i=0
    for L in Ls:
        t_vals_15 = np.zeros(4)
        t_sigma_15 = np.zeros(4)
        t_vals_23 = np.zeros(4)
        t_sigma_23 = np.zeros(4)
        for nproc in nps:
            filename=f"timeMPI_L{L}_nproc{nproc}.txt"
            fil = open(filename, 'r')
            lines = fil.readlines()
            T15 = []
            T23 = []
            for line in lines:
                split = line.split('=')
                t = split[1].split()[0]
                temp = split[2].split()[0].strip()
                if temp == '1.5':
                    T15.append(float(t))
                elif temp == '2.3':
                    T23.append(float(t))
            fil.close()
            t_vals_15[nproc-1] = np.mean(T15)
            t_sigma_15[nproc-1] = np.std(T15)/np.sqrt(len(T15))
            t_vals_23[nproc-1] = np.mean(T23)
            t_sigma_23[nproc-1] = np.std(T23)/np.sqrt(len(T23))
        axs[i,0].errorbar(nps,t_vals_15,yerr=t_sigma_15,fmt='xb')
        axs[i,0].plot(np.linspace(1,4,1000),t_vals_15[0]/np.linspace(1,4,1000),'b--')
        
        axs[i,1].errorbar(nps,t_vals_23,yerr=t_sigma_23,fmt='xb')
        axs[i,1].plot(np.linspace(1,4,1000),t_vals_23[0]/np.linspace(1,4,1000),'b--')
        i+=1
    
    for ax in axs.flat:
        ax.set(xlabel='Number of cores', ylabel='Mean time [s]')
    for ax in axs.flat:
        ax.label_outer()
    axs[0, 0].set_title('T = 1.5, L = 40')
    axs[0, 1].set_title('T = 2.3, L = 40')
    axs[1, 0].set_title('T = 1.5, L = 100')
    axs[1, 1].set_title('T = 2.3, L = 100')
    plt.tight_layout()
    plt.savefig("MPItime.pdf")
    plt.show()


plot_mpi()
#plot_mpi_timing()
"""
Plotting of the energy and magnetization over time

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
"""