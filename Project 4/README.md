# Project 4  
Numerical solutions of the two dimensional Ising model

## Program descriptions  
- **isingmodel.cpp**  
Contains the functions used for the Metropolis algorithm, and initialization of the system

- **main.cpp**  
File used to run the non parallelized simulations  

- **mainMPI.cpp**  
File used to run the parallelized simulations for larger system sizes  

- **makefile**  
Compiles all programs

- **run_mpi.sh**  
Bash script used to run the simulations for L=40 to L=100.  
(Used on a computer with 4 cores, and Ubuntu 18.04)  

- **time_mpi.sh**  
Bash script used to time the parallelized version of the program with different amounts of cores.  
Used on a comuter with 4 cores.  
