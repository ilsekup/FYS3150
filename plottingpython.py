import matplotlib.pyplot as plt
import numpy as np
import sys, os

filename = sys.argv[1] #filename in terminal when running program
fil=open(filename, 'r') #reading file
line1 = fil.readline()
line2 = fil.readline()
max_error = float(line1.split("=")[1])
T = float(line2.split("=")[1])
fil.readline()
fil.readline()
lines=fil.readlines() #splitting in lines
N = len(lines)
v=np.zeros(N)#[]  #making separate lists for the different collums
u=np.zeros(N)#[]
x=np.zeros(N)#[]
for num,i in enumerate(lines): #putting everything in the right list
    text = i.split(",")
    x[num] = float(text[0])
    u[num] = float(text[1])
    v[num] = float(text[2])
fil.close()

x = np.array(x)
u = np.array(u)
v = np.array(v)

fil2=open("run_info_standard.txt", 'r')
lines = fil2.readlines()
num = len(lines)
N = []; eps = []; time = []
number_of_runs = []
for i in range(num):
    text = lines[i].split("=")
    N_ = int(text[1].split()[0])
    time_ = float(text[2].split()[0])
    eps_ = float(text[3].split()[0])
    if N_ in N:
        index = N.index(N_)
        number_of_runs[index]+=1
        eps[index] += eps_
        time[index] += time_
    else:
        N.append(N_)
        eps.append(eps_)
        time.append(time_)
        number_of_runs.append(1)

fil2.close()
t_standard = np.array(time)
err_standard = np.array(eps)
N_standard = np.array(N)
N_runs_standard = np.array(number_of_runs)
t_standard /= N_runs_standard; err_standard /= N_runs_standard

fil3=open("run_info_specialized.txt", 'r')
lines = fil3.readlines()
num = len(lines)
print("num = ", num)
N = []; eps = []; time = []
number_of_runs = []
for i in range(num):
    text = lines[i].split("=")
    N_ = int(text[1].split()[0])
    time_ = float(text[2].split()[0])
    eps_ = float(text[3].split()[0])
    if N_ in N:
        index = N.index(N_)
        number_of_runs[index]+=1
        eps[index] += eps_
        time[index] += time_
    else:
        N.append(N_)
        eps.append(eps_)
        time.append(time_)
        number_of_runs.append(1)

fil3.close()

t_specialized = np.array(time)
err_specialized = np.array(eps)
N_specialized = np.array(N)
N_runs_specialized = np.array(number_of_runs)
print(N_runs_specialized)
t_specialized /= N_runs_specialized; err_specialized /= N_runs_specialized

fil4=open("run_info_LU.txt", 'r')
lines = fil4.readlines()
num = len(lines)
print("num = ", num)
N = []; eps = []; time = []
number_of_runs = []
for i in range(num):
    text = lines[i].split("=")
    N_ = int(text[1].split()[0])
    time_ = float(text[2].split()[0])
    eps_ = float(text[3].split()[0])
    if N_ in N:
        index = N.index(N_)
        number_of_runs[index]+=1
        eps[index] += eps_
        time[index] += time_
    else:
        N.append(N_)
        eps.append(eps_)
        time.append(time_)
        number_of_runs.append(1)

fil4.close()

t_LU = np.array(time)
err_LU = np.array(eps)
N_LU = np.array(N)
N_runs_LU = np.array(number_of_runs)
print(N_runs_LU)
t_LU /= N_runs_LU; err_LU /= N_runs_LU

plt.figure(figsize = (7,5))
plt.title("Time used")
plt.loglog(N_standard,t_standard, label = "Standard Thomas")
plt.loglog(N_specialized,t_specialized, label = "Specialized Thomas")
plt.legend(fontsize = 12)
plt.xlabel("N",fontsize = 14)
plt.ylabel("Time [s]",fontsize = 14)
plt.grid()
plt.savefig("elapsed_time.png")

plt.figure(figsize = (7,5))
plt.title("Error in numerical simulation",fontsize = 16)
plt.semilogx(N_standard,err_standard, label = "Standard Thomas")
plt.semilogx(N_specialized,err_specialized, label = "Specialized Thomas")
plt.legend(fontsize = 12)
plt.xlabel("N",fontsize = 14)
plt.ylabel(r"$\log_{10}$ of max error",fontsize = 14)
plt.grid()
plt.savefig("error.png")

plt.figure()

def exact(x):
    return 1 - (1- np.exp(-10))*x - np.exp(-10*x)

xi = np.linspace(0, 1, 100)
plt.plot(x, v, label="calculated")
plt.plot(xi, exact(xi), 'r--', label="exact")
plt.xlabel('x values')
plt.ylabel('y values')
plt.title('Result of the numerical integration')
plt.grid()
plt.legend(fontsize = 12)

plt.figure(figsize = (7,5))
print(N_LU)
plt.title("Time used LU method",fontsize = 16)
plt.loglog(N_LU,t_LU)
plt.xlabel("N",fontsize = 14)
plt.ylabel("Time [s]",fontsize = 14)
plt.grid()
plt.savefig("elapsed_time_LU.png")

plt.show()