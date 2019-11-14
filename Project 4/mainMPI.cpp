#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "isingmodel.h"
#include "isingmodel.cpp"
#include"mpi.h"

using namespace std;
int main(int argc, char* argv[]){
  char *outfilename;
  int n, my_rank, numprocs;
  double mc;
  double w[17], average[5], total_average[5], initial_T, final_T, E, M, T_step;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if(my_rank==0 && argc > 1){
    outfilename = argv[1];
    ofile.open(outfilename);
  }
  n = 20;
  mc = 100000;
  int no_intervalls = mc/numprocs;
  int myloop_begin = my_rank*no_intervalls + 1;
  int myloop_end = (my_rank +1)*no_intervalls;
  if((my_rank == numprocs -1) && (myloop_end < mc)){
    myloop_end = mc;
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&initial_T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&final_T,1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&T_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  int **spin_matrix = (int**) matrix(n, n, sizeof(int));
  long startpoint = -1-my_rank;
  for(double temp = initial_T; temp <= final_T; temp += T_step) {
    E = M = 0;
    for(int de=-8; de <= 8; de++) w[de+8] = 0;
    for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/temp);
    for (int i = 0; i < 5; i++) average[i] = 0.;
    for(int i = 0; i < 5; i++) total_average[i] = 0.;
    initialize(n,  temp, spin_matrix, E, M, false);  //false if we want random spins, true for ordered
    for(int cycles=1; cycles <= mc; cycles++){
      metropolis(n, startpoint, spin_matrix, E, M, w);
      average[0] += E;
      average[1] += E*E;
      average[2] += M;
      average[3] += M*M;
      average[4] += fabs(M);
     }
     for(int i = 0; i <5; i++){
       MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     }
     if(my_rank==0){
       writingfunc(n, mc, temp, total_average);
     }
  }
  free_matrix((void **) spin_matrix);
  ofile.close();
  MPI_Finalize();
  return 0;
}