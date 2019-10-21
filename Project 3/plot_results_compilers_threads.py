
import numpy as np
import matplotlib.pyplot as plt

exact_answer = 5*np.pi**2/256
#plot_results.py and plot_results_threadsvstime rewritten to fit compiler flags
#+ threads
filenames = ["estimates_MonteCarlo1.txt","estimates_MonteCarlo4.txt", "estimates_MonteCarloO2.txt", "estimates_MonteCarloO3.txt"]
#Change the name of the files and labels to fit your plotting needs
labels1 = ["Normal", "4 threads", "4 threads + -O2", "4 threads + -O3" ]
labels2 = ["Normal (IS)", "4 threads (IS)", "4 threads + -O2 (IS)", "4 threads + -O3 (IS)" ]
threads = np.array([1,2,3,4]) #number of threads used
plt.figure(figsize = (7,5))

for k in range(len(filenames)):
    fil = open(filenames[k], 'r') #reading file

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
        sigma_avg[i] = np.sum(sigma[bool_arr])/count[i]
        mean_error_MC[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count[i]
        t_avg[i] = np.sum(t[bool_arr])/count[i]
        act_sig[i] = np.sqrt(np.sum((exact_answer - I[bool_arr])**2)/count[i])
    plt.plot(np.log10(t_avg),np.log10(mean_error_MC)  ,"o",label = labels1[k] ,markersize = 6)


# Repeat for importance sampling

filenames = ["estimates_MonteCarlo_imps1.txt","estimates_MonteCarlo_imps4.txt", "estimates_MonteCarlo_impsO2.txt", "estimates_MonteCarlo_impsO3.txt"]

for k in range(len(filenames)):
    fil = open(filenames[k], 'r') #reading file

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
    act_sig_imps = np.zeros(len(Ns_imps))

    for i,n in enumerate(Ns_imps):
        bool_arr = N == n
        sigma_imps_avg[i] = np.sum(sigma[bool_arr])/count_imps[i]
        mean_error_MC_imps[i] = np.sum(np.abs(I[bool_arr]-exact_answer))/count_imps[i]
        t_imps_avg[i] = np.sum(t[bool_arr])/count_imps[i]
        act_sig_imps[i] = np.sqrt(np.sum((exact_answer - I[bool_arr])**2)/count_imps[i])

    plt.plot(np.log10(t_imps_avg),np.log10(mean_error_MC_imps) ,"o",label = labels2[k] ,markersize = 6)


plt.xlabel(r"$\log_{10}$ time [s]", fontsize = 14)
plt.ylabel(r"$\log_{10}$ mean error", fontsize = 14)
plt.title("Mean error as a function of mean time",fontsize = 16)
plt.grid()
plt.legend(fontsize = 12)
plt.savefig("Parallelization + compilation.png")
plt.show()
