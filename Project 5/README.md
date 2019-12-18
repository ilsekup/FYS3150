# Project 5  
Numerical partial differential equation solver 

## Program descriptions  
- **PDEsolver.cpp**  
Contains the functions used for the explicit, implicit and crank-nicolson scheme in 1D and the implicit and  explicit in 2D in addition to the geological case

- **main.cpp**  
File used to run the functions, needs to be loaded with armadillo and PDEsolver.cpp
 

- **makefile**  
Compiles all programs and runs testfunctions 

- **plottingdiff.py**  
Makes animation of one dimensional data.

- **2dplotting.py**  
Makes contour plots of the 2D data. Needs ffmpeg installed in order to function properly and make an mp4 file, without this installed the frames will not be removed and a temp folder with a lot of files will be stored on computer. 

