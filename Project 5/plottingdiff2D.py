

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


filename2D = "explicit.txt"
u_t = []

with open(filename2D, 'r') as infile:
    infile.readline()
    u_yx_mat = []
    for line in infile:
        u_yx_vec = []
        for i in range(len(line.split())):
            u_yx_vec.append(float(line.split()[i]))
        u_yx_mat.append(u_yx_vec)
    else:
        infile.readline()



"""
n = 10
x = np.linspace(0,1,n+1)
y = np.linspace(0,1,n+1)

X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize = (16,8))
ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
line, = ax.plot([], [], lw=2)

plt.xlabel("x-position")
plt.ylabel("y-position")
plt.title("Temporal evolution of density u as a function of the x and y-position")
plt.show()
"""
