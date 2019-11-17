#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
#include "isingmodel.h"
#include "mpi.h"
#include "lib.h"

using namespace std;
int main(int argc, char* argv[]){
  using namespace std::chrono;
  char *outfilename;
  int n, my_rank, numprocs;
  int mc;
  double w[17], average[5], total_average[5], initial_T, final_T, E, M, T_step;
  initial_T = 2.0;
  final_T = 2.3;
  T_step = 0.01;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  //if(my_rank==0 && argc > 1){
    outfilename = argv[1];
    ofstream ofile;
    ofile.open(outfilename);
  //}
  mc = atof(argv[2]);
  n = atof(argv[3]);
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
  long seed = time(NULL) + my_rank * 1000;
  int mc_per_proc = (int) mc/numprocs;
  high_resolution_clock::time_point t1, t2;
  duration<double, ratio<1,1>> t;

  if(my_rank == 0){
    mc_per_proc += mc%numprocs;
  }

  for(double temp = initial_T; temp <= final_T; temp += T_step) {
    if(my_rank == 0){
      t1 = high_resolution_clock::now();
    }
    E = M = 0;
    for(int de=-8; de <= 8; de++) w[de+8] = 0;
    for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/temp);
    for (int i = 0; i < 5; i++) average[i] = 0.;
    for(int i = 0; i < 5; i++) total_average[i] = 0.;
    initialize(n,  temp, spin_matrix, E, M, seed, false);  //false if we want random spins, true for ordered
    for(int cycles=1; cycles <= mc_per_proc; cycles++){
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
     //save results and print time and current progress
     if(my_rank==0){
       writingfunc(n, mc, temp, total_average, ofile);
       t2 = high_resolution_clock::now();
       t = t2-t1;
       cout << "Time spent on MC = " << t.count() << " seconds, T = " << temp << " L = " << n <<endl;
     }
  }
  free_matrix((void **) spin_matrix);
  ofile.close();
  MPI_Finalize();
  return 0;
}