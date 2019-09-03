import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1] #filename in terminal when running program 
fil=open(filename, 'r') #reading file
lines=fil.readlines() #splitting in lines
N = len(lines)
v=np.zeros(N)#[]  #making separate lists for the different collums
u=np.zeros(N)#[]
x=np.zeros(N)#[]
for num,i in enumerate(lines): #putting everything in the right list
    text = i.split()
    x[num] = float(text[1])
    u[num] = float(text[3])
    v[num] = float(text[-1])
fil.close()

x = np.array(x)
u = np.array(u)
v = np.array(v)

def exact(x):
    return 1 - (1- np.exp(-10))*x - np.exp(-10*x)

xi = np.linspace(0, 1, 100)
#print(exact(xi))
#print(u)
#plt.plot(x, u)
plt.plot(x, v, label="calculated") #neither u or v look like they should, at all
plt.plot(xi, exact(xi),label="exact")
plt.xlabel('x values')
plt.ylabel('y values')
plt.legend()
plt.show()
