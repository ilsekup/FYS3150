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

def exact(x):
    return 1 - (1- np.exp(-10))*x - np.exp(-10*x)

xi = np.linspace(0, 1, 100)
plt.plot(x, v, label="calculated")
plt.plot(xi, exact(xi), 'r--', label="exact")
plt.xlabel('x values')
plt.ylabel('y values')
plt.title('Result of the numerical integration')
plt.grid()
plt.legend()
plt.show()
