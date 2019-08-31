import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1] #filename in terminal when running program 
fil=open(filename, 'r') #reading file
lines=fil.readlines() #splitting in lines
v=[]  #making separate lists for the different collums
u=[]
x =[]
for i in lines: #putting everything in the right list
    x.append(i.split()[1])
    u.append(i.split()[3])
    v.append(i.split()[-1])
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
plt.plot(x, v) #neither u or v look like they should, at all
plt.plot(xi, exact(xi))
plt.xlabel('x values')
plt.ylabel('y values')
plt.show()
