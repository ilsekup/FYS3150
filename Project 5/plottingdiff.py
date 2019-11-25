import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


filename1 = "explicit.txt"
file1 = open(filename1, 'r')
lines1 = file1.readlines()

u_xx = [] # nested list containing u for each x at every t-iteration
for line in lines1:    #1 2 3 --- t lines
    u_xx_temp = []
    for i in range(len(line.split())): # append u for each x and time step
        u_xx_temp.append(float(line.split()[i]))
    u_xx.append(u_xx_temp) # append the total u for each x for a given t


t = 0
x = np.linspace(0,1, len(u_xx[0]))

#Making an animation to see temporal evolution
#when running the file you get from PDEsolver.cpp
#Try n = 100, t = 20 000 and running this code
fig = plt.figure(figsize = (16,8))
ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(0,1, len(u_xx[0]))
    y = u_xx[i]
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2000, interval=200, blit=True)
plt.xlabel("x-position")
plt.ylabel("u density")
plt.title("Temporal evultion of density u as a function of the x-position")
plt.show()
