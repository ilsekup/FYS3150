import numpy as np
import matplotlib.pyplot as plt
import sys

exact_answer = 5*np.pi**2/256

filename = "estimates_MonteCarlo.txt"
fil = open(filename, 'r') #reading file

lines = fil.readlines() #splitting in lines
N_data = len(lines)

# Initialize lists
N = np.zeros(N_data, dtype = int)
I = np.zeros(N_data)
sigma = np.zeros(N_data)
t = np.zeros(N_data)

for i,line in enumerate(lines): #Putting the eigenvectors into arrays
    text = line.split('=')
    _N = int(text[1].split()[0])
    _I = float(text[2].split()[0])
    _sigma = float(text[3].split()[0])
    _t = float(text[4].split()[0])

    N[i] = _N; I[i] = _I; sigma[i] = _sigma; t[i] = _t
    

fil.close()

Ns, count = np.unique(N,False,False,True)
I_avg = np.zeros(len(Ns))
sigma_avg = np.zeros(len(Ns))
t_avg = np.zeros(len(Ns))
mean_error_MC = np.zeros(len(Ns))

for i,n in enumerate(Ns):
    bool_arr = N == n
    print(f"num {n} normal  =  "  + str(count[i]))
    sigma_avg[i] = np.sum(sigma[bool_arr])/count[i]
    mean_error_MC[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count[i]
    t_avg[i] = np.sum(t[bool_arr])/count[i]

# Repeat for importance sampling
filename = "estimates_MonteCarlo_imps.txt"
fil = open(filename, 'r') #reading file

lines = fil.readlines() #splitting in lines
N_data = len(lines)

# Initialize lists
N = np.zeros(N_data, dtype = int)
I = np.zeros(N_data)
sigma = np.zeros(N_data)
t = np.zeros(N_data)

for i,line in enumerate(lines): #Putting the eigenvectors into arrays
    text = line.split('=')
    _N = int(text[1].split()[0])
    _I = float(text[2].split()[0])
    _sigma = float(text[3].split()[0])
    _t = float(text[4].split()[0])

    N[i] = _N; I[i] = _I; sigma[i] = _sigma; t[i] = _t
    

fil.close()

Ns_imps, count_imps = np.unique(N,False,False,True)
sigma_imps_avg = np.zeros(len(Ns_imps))
mean_error_MC_imps = np.zeros(len(Ns_imps))
t_imps_avg = np.zeros(len(Ns_imps))

for i,n in enumerate(Ns_imps):
    bool_arr = N == n
    print(f"num {n} imps  =  "  + str(count_imps[i]))
    sigma_imps_avg[i] = np.sum(sigma[bool_arr])/count_imps[i]
    mean_error_MC_imps[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count_imps[i]
    t_imps_avg[i] = np.sum(t[bool_arr])/count_imps[i] 

# Find gauss quadrature values

filename = "estimates_gauss_legendre.txt"
fil = open(filename, 'r') #reading file

lines = fil.readlines() #splitting in lines
N_data = len(lines)

# Initialize lists
N = np.zeros(N_data, dtype = int)
I = np.zeros(N_data)
t = np.zeros(N_data)

for i,line in enumerate(lines): #Putting the eigenvectors into arrays
    text = line.split('=')
    _N = int(text[1].split()[0])
    _I = float(text[2].split()[0])
    _t = float(text[3].split()[0])
    N[i] = _N; I[i] = _I; t[i] = _t
    

fil.close()

Ns_leg, count_leg = np.unique(N,False,False,True)
mean_error_leg = np.zeros(len(Ns_leg))
t_leg_avg = np.zeros(len(Ns_leg))

for i,n in enumerate(Ns_leg):
    bool_arr = N == n
    print(f"num {n} leg  =  "  + str(count_leg[i]))
    mean_error_leg[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count_leg[i]
    t_leg_avg[i] = np.sum(t[bool_arr])/count_leg[i] 

filename = "estimates_gauss_laguerre.txt"
fil = open(filename, 'r') #reading file

lines = fil.readlines() #splitting in lines
N_data = len(lines)

# Initialize lists
N = np.zeros(N_data, dtype = int)
I = np.zeros(N_data)
t = np.zeros(N_data)

for i,line in enumerate(lines): #Putting the eigenvectors into arrays
    text = line.split('=')
    _N = int(text[1].split()[0])
    _I = float(text[2].split()[0])
    _t = float(text[3].split()[0])

    N[i] = _N; I[i] = _I; t[i] = _t
    

fil.close()

Ns_lag, count_lag = np.unique(N,False,False,True)
mean_error_lag = np.zeros(len(Ns_lag))
t_lag_avg = np.zeros(len(Ns_lag))

for i,n in enumerate(Ns_lag):
    bool_arr = N == n
    print(f"num {n} leg  =  "  + str(count_leg[i]))
    mean_error_lag[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count_lag[i]
    t_lag_avg[i] = np.sum(t[bool_arr])/count_lag[i]

plt.figure(figsize = (7,5))
plt.title("The mean error, as function of integration points",fontsize = 16)
plt.xlabel(r"$\log_{10}$ N", fontsize = 14)
plt.ylabel(r"$\log_{10}$ mean error", fontsize = 14)
plt.plot(np.log10(Ns),np.log10(mean_error_MC),label = "Brute force MC")
plt.plot(np.log10(Ns_imps),np.log10(mean_error_MC_imps),label = "Importance sampling MC")


plt.legend(fontsize = 12)
plt.savefig("errors_MC.png")
plt.show()

plt.figure(figsize = (7,5))
plt.xlabel(r"$\log_{10}$ N", fontsize = 14)
plt.ylabel(r"$\log_{10}$ time [s]", fontsize = 14)
plt.title("The time used, as function of integration points",fontsize = 16)
plt.plot(np.log10(Ns),np.log10(t_avg),label = "Brute force")
plt.plot(np.log10(Ns_imps),np.log10(t_imps_avg),label = "Importance sampling")

plt.legend(fontsize = 12)
plt.savefig("time_MC.png")

plt.figure(figsize = (7,5))
plt.xlabel(r"$\log_{10}$ time [s]", fontsize = 14)
plt.ylabel(r"$\log_{10}$ mean error", fontsize = 14)
plt.title("Mean error, as function of time",fontsize = 16)
plt.plot(np.log10(t_avg),np.log10(mean_error_MC), "x",label = "Brute force MC")
plt.plot(np.log10(t_imps_avg),np.log10(mean_error_MC_imps), "x",label = "Importance sampling MC")
plt.plot(np.log10(t_leg_avg),np.log10(mean_error_leg), "x",label = "Legendre MC")
plt.plot(np.log10(t_lag_avg),np.log10(mean_error_lag), "x",label = "Laguerre MC")

plt.legend(fontsize = 12)
plt.savefig("errorVStime.png")

plt.show()