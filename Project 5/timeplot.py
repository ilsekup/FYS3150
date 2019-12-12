import numpy as np
import matplotlib.pyplot as plt

filename = "Time.txt"

file = open(filename,'r')
lines = file.readlines()

n = []
thread1 = []
thread2 = []
thread3 = []
thread4 = []

for line in lines:
    n.append(float(line.split()[0]))
    thread1.append(float(line.split()[1]))
    thread2.append(float(line.split()[2]))
    thread3.append(float(line.split()[3]))
    thread4.append(float(line.split()[4]))


#plotting

plt.figure(figsize=(8,6))
plt.title("Time spent vs size of matrices in 2D explicit scheme, # t steps = 3 n^2",fontsize = 14)
plt.plot(n,thread1, label=" 1 thread")
plt.plot(n,thread2, label=" 2 threads")
plt.plot(n,thread3, label=" 3 threads")
plt.plot(n,thread4, label=" 4 threads")
plt.xlabel("n, size of matrix (n x n)", fontsize = 14)
plt.ylabel("Time spent on code for explicit 2D [s]", fontsize = 14)
plt.legend()
plt.tight_layout()
plt.show()
