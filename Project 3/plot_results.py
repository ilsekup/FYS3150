import numpy as np
import matplotlib.pyplot as plt
from scipy import stats 

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

for i,line in enumerate(lines):
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
act_sig = np.zeros(len(Ns))

for i,n in enumerate(Ns):
    bool_arr = N == n
    print(f"num {n} normal  =  "  + str(count[i]))
    sigma_avg[i] = np.sum(sigma[bool_arr])/count[i]
    mean_error_MC[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count[i]
    t_avg[i] = np.sum(t[bool_arr])/count[i]
    act_sig[i] = np.sqrt(np.sum((exact_answer - I[bool_arr])**2)/count[i])


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

for i,line in enumerate(lines):
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
act_sig_imps = np.zeros(len(Ns_imps))

for i,n in enumerate(Ns_imps):
    bool_arr = N == n
    print(f"num {n} imps  =  "  + str(count_imps[i]))
    sigma_imps_avg[i] = np.sum(sigma[bool_arr])/count_imps[i]
    mean_error_MC_imps[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count_imps[i]
    t_imps_avg[i] = np.sum(t[bool_arr])/count_imps[i] 
    act_sig_imps[i] = np.sqrt(np.sum((exact_answer - I[bool_arr])**2)/count_imps[i])

# Find gauss quadrature values

filename = "estimates_gauss_legendre.txt"
fil = open(filename, 'r') #reading file

lines = fil.readlines() #splitting in lines
N_data = len(lines)

# Initialize lists
N = np.zeros(N_data, dtype = int)
I = np.zeros(N_data)
t = np.zeros(N_data)

for i,line in enumerate(lines):
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

for i,line in enumerate(lines):
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

slope, intercept, _,_, std_err = stats.linregress(np.log10(Ns),np.log10(sigma_avg))
slope_imps, intercept_imps, _,_, std_err_imps = stats.linregress(np.log10(Ns),np.log10(sigma_imps_avg))

coeffs_bruteforce = np.polyfit(np.log10(Ns),np.log10(sigma_avg),1)
coeffs_imps = np.polyfit(np.log10(Ns),np.log10(sigma_imps_avg),1)

print(f"slope for brute force = {slope:.2f} ± {std_err:.2f}")
print(f"slope for importance sampling = {slope_imps:.3f} ± {std_err_imps:.3f}")


plt.figure(figsize = (7,5))
plt.title(r"The mean $\sigma$ plotted against integration points",fontsize = 16)
plt.xlabel(r"$\log_{10}$ N", fontsize = 14)
plt.ylabel(r"$\log_{10} \, \overline{\sigma}$", fontsize = 14)

plt.plot(np.log10(Ns),np.log10(sigma_avg),'bx',label = r"Brute force MC ($\infty = 2$)",markersize = 8)
plt.plot(np.log10(Ns_imps),np.log10(sigma_imps_avg),'rx',label = "Importance sampling MC",markersize = 8)

plt.plot(np.log10(Ns),np.log10(act_sig),'bo')
plt.plot(np.log10(Ns_imps),np.log10(act_sig_imps),'ro')

plt.plot(np.log10(Ns),np.polyval(coeffs_bruteforce,np.log10(Ns)), 'b--')
plt.plot(np.log10(Ns),np.polyval(coeffs_imps,np.log10(Ns)), 'r--')

plt.grid()
plt.legend(fontsize = 12)
plt.savefig("errors_MC.png")

plt.figure(figsize = (7,5))
plt.xlabel(r"$\log_{10}$ N", fontsize = 14)
plt.ylabel(r"$\log_{10}$ time [s]", fontsize = 14)
plt.title("The time used plotted against of integration points",fontsize = 16)
plt.plot(np.log10(Ns),np.log10(t_avg),'x',label = r"Brute force ($\infty = 2$)",markersize = 8)
plt.plot(np.log10(Ns_imps),np.log10(t_imps_avg),'x',label = "Importance sampling",markersize = 8)

plt.legend(fontsize = 12)
plt.grid()

plt.savefig("time_MC.png")

plt.figure(figsize = (7,5))
plt.xlabel(r"$\log_{10}$ time [s]", fontsize = 14)
plt.ylabel(r"$\log_{10}$ mean error", fontsize = 14)
plt.title("Mean error, as function of time",fontsize = 16)
plt.plot(np.log10(t_avg),np.log10(mean_error_MC), "x",label = r"Brute force MC ($\infty = 2$)",markersize = 8)
plt.plot(np.log10(t_imps_avg),np.log10(mean_error_MC_imps), "x",label = "Importance sampling MC",markersize = 8)
plt.plot(np.log10(t_leg_avg),np.log10(mean_error_leg), "x",label = r"Legendre GQ ($\infty = 2$)",markersize = 8)
plt.plot(np.log10(t_lag_avg),np.log10(mean_error_lag), "x",label = "Laguerre GQ",markersize = 8)

plt.grid()
plt.legend(fontsize = 12)
plt.savefig("errorVStime.png")

plt.show()