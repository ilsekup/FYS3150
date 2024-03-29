import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys, os, shutil
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
This script either makes a gif for the non specific case, or the geological case
as sys.argv[1] give the filename of the file to plot
For the non-specific case, set plotfunc (sys.argv[2]) to plot_square
For the geological case, set plotfunc to plot_geo
"""

def plot_square(X,Y,vals,N,time,i,max_val):
    """
    Plots contour plot for the square case (not representing any physical case)
    """
    levels = 128
    fig = plt.figure(figsize = (10,7))
    plt.contourf(X,Y,vals,levels)
    plt.xlabel("x", fontsize = 14)
    plt.ylabel("y", fontsize = 14)
    plt.title(f"2D diffusion, n = {N-2}, t = {time:.4f}",fontsize=16)
    cbar = plt.colorbar(ticks=np.linspace(0,1,11))
    cbar.ax.set_ylabel('Temperature (reduced units)', rotation=270, labelpad=15)
    plt.tight_layout()
    plt.savefig(f"frame_{i:04d}.png")
    plt.close(fig)

def plot_geo(X,Y,vals,nx,time,i,max_val):
    """
    Plots contour plot for the geological case
    """
    levels = np.linspace(0,max_val*1573.0,128)

    # Rescale units
    X_plot = (1-X)*120.0
    Y_plot = Y*120.0
    vals_plot = vals*1573.0
    vals_plot -= 273.15
    time_use = time*0.6392
    base_vals = np.array([np.linspace(0.1787,1,nx)[::-1] for i in range(np.max(X.shape))]).T
    base_vals *= 1573.0
    base_vals -= 273.15
    vals_plot -= base_vals
    fig = plt.figure(figsize = (10,7))
    plt.contourf(Y_plot,X_plot,vals_plot.T,levels,extend='both')
    plt.axis("equal")
    plt.xlabel("Width [km]", fontsize = 14)
    plt.ylabel("Depth [km]", fontsize = 14)
    plt.title(rf"Heatflow in lithosphere, $n_x$ = {nx-2}, t = {0.6392*time:.4f} Gyr",fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Temperature difference [${}^\circ$ C]', rotation=270, labelpad=15)
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.savefig(f"frame_{i:04d}.png")
    plt.close(fig)

def make_gif(filename,plotting_func):
    """
    Plots the first first tasks (with square matrices) and stores them in temp
    """
    anim_time = 7 #s
    fps = 25
    tot_frames = int(anim_time*fps)

    infofile = open("runinfo.txt", 'r')
    line = infofile.readline().split('=')
    dt = float(line[1].split()[0].strip())
    nx = int(line[2].split()[0].strip())+2
    ny = int(line[3].split()[0].strip())+2
    t_stop = float(line[4].split()[0].strip())
    infofile.close()

    dims = [nx, ny]
    arr = np.zeros(dims)

    temppath = 'temp'

    if not os.path.exists(temppath):
        os.mkdir(temppath)

    with open(filename, 'r') as file:
        next(file)

        i = 0
        num_frames = 0
        for line in file:
            if line == '\n':
                np.save(f"temp/frame{num_frames:05d}.npy",arr)
                num_frames += 1
                print(f"Timestep {num_frames} saved")
                arr = np.zeros(dims)
                
                i = 0
            else:
                arr[i] = np.fromstring(line, sep = ' ')
                i+=1

    n_skip = int(num_frames/tot_frames)
    use_frames = np.arange(0,num_frames)[::n_skip]

    frames = []

    for i,frame in enumerate(use_frames):
        print(f"Loading frame {i} of {len(use_frames)}")
        arr = np.load(f"temp/frame{frame:05d}.npy")
        frames.append(arr.copy())

    # For geo-plotting
    base_vals = np.array([np.linspace(0.1787,1,nx)[::-1] for i in range(ny)]).T
    max_val = 0
    for frame in frames:
        if np.max(frame-base_vals) > max_val:
            max_val = np.max(frame-base_vals)

    x = np.linspace(0,1,nx)
    ymax = 1*ny/nx
    y = np.linspace(0,ymax,ny)
    X,Y=np.meshgrid(x,y)
    levels = 128


    os.chdir(temppath)

    for i,frame in enumerate(frames):
        time = n_skip*i*dt;
        print(f"Plotting frame {i} of {len(frames)}")
        plotting_func(X,Y,frame,nx,time,i,max_val)

    file_head = filename.split(".")[0]
    os.system(f"ffmpeg -r {fps} -i frame_%04d.png -vcodec libx264 -crf 20 {file_head}.mp4")
    os.chdir("..")
    shutil.move(f"temp/{file_head}.mp4", f"Results/{file_head}.mp4")
    transfer_frames = [0,20,80,i]
    for k,j in enumerate(transfer_frames):
        shutil.move(f"temp/frame_{j:04d}.png", f"Results/{file_head}_{k}.png")

    shutil.rmtree('temp')

filename = sys.argv[1]
plotfunc_string = sys.argv[2]

if plotfunc_string=="plot_square":
    plt_func = plot_square
elif plotfunc_string=="plot_geo":
    plt_func = plot_geo
else:
    print("Second command line argument must be either 'plot_square' or 'plot_geo'")
    sys.exit(1)

make_gif(filename,plt_func)