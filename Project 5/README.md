# Project 5  
Numerical partial differential equation solver 

## Program descriptions  
- **PDEsolver.cpp**  
Contains the functions used for the explicit, implicit and crank-nicolson scheme in 1D and the implicit and  explicit in 2D in addition to the geological case

- **main.cpp**  
File used to run the functions, needs to be loaded with armadillo and PDEsolver.cpp

- **Geomain.cpp**  
File used to simulate the three geological cases, must be loaded with armadillo.
In the report this was run with n=100, and it automatically runs the last case to 3 Gyr (as time evolution is more important for this case)
Note that these parameters produce a 22 GB textfile. The value for N could therefore be lowered to save space.  

- **2dplotting.py**  
This plots filled countour plots of the two dimensional simulations and saves an mp4 animation in the Results folder. Must be run with two command line arguments, where the first is the filename, and the second is the plotting function to be used. This must be either "plot_square" for the first cases (with square matrices and reduced units) or "plot_geo" if the file represents a simulation of the geological case. Note that this program saves copies of each timestep to avoid storing them in memory. So if the textfile from the simulation is large, this program will temporarily use much storage space. This is however deleted at the end of the program. It also requires ffmpeg for the animation.
 
- **makefile**  
Compiles all programs and runs testfunctions 

- **plottingdiff.py**  
Makes animation of one dimensional data.

- **2dplotting.py**  
Makes contour plots of the 2D data. Needs ffmpeg installed in order to function properly and make an mp4 file, without this installed the frames will not be removed and a temp folder with a lot of files will be stored on computer. 

