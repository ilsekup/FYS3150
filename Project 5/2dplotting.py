import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys, os
from mpl_toolkits.axes_grid1 import make_axes_locatable

filename = sys.argv[1]
file = open(filename, 'r')

line1 = file.readline()
N = int(line1.strip())
file.readline()
lines = file.readlines()


dims = [N, N]

arr = np.zeros(dims)

values = []

i = 0
for line in lines:
    if line == '\n':
        print("yes")
        arr = np.zeros(dims)
        values.append(arr)
        i = 0
    else:
        arr[i] = np.fromstring(line, sep = ' ')
        i+=1

xmin, xmax, ymin, ymax = 0, 1, 0, 1
x = y = np.linspace(0,1,N)
X,Y=np.meshgrid(x,y)
levels = 60

anim_time = 10 #s
fps = 60
tot_frames = anim_time*fps

n_skip = int(round(len(values)/tot_frames))

frames = values[::n_skip]

temppath = 'temp'

if not os.path.exists(temppath):
    os.mkdir(temppath)

os.chdir(temppath)

for i,frame in enumerate(frames):
    fig = plt.figure(figsize = (10,7))
    plt.contourf(X,Y,frame,levels)
    plt.xlabel("x/L", fontsize = 14)
    plt.ylabel("y/L", fontsize = 14)
    plt.colorbar()
    plt.savefig(f"frame_{i:04d}.png")
    print(i)
    plt.close(fig)

os.system(f"convert -delay {100/fps} -loop 0 *.png output.gif")
"""
fig = plt.figure(figsize = (10,7))
ax = plt.axes(xlim=(xmin,xmax),ylim = (ymin,ymax))
def init():
    plane = ax.contourf(X,Y,frames[0],levels)
    plt.axis("equal")
    return plane

def anim(i):
    plane = ax.contourf(X,Y,frames[i],levels)
    plt.axis("equal")
    print(i)
    return plane

ani = animation.FuncAnimation(fig,anim,frames = int(fps*anim_time),init_func = init, interval = 1/fps, repeat = False)
writer = animation.FFMpegWriter()
ani.save("animation.mp4")
"""