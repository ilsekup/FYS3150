#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "isingmodel.h"
#include "isingmodel.cpp"

using namespace std;
int main(int argc, char* argv[]){
  char *outfilename;
  int n, mc;
  double w[17], average[5], initial_T, final_T, E, M, T_step;
  outfilename = argv[1];
  ofile.open(outfilename);
  cout << "Number of spins: " << endl;
  cin >> n;
  cout << "Number of monte carlo cycles" << endl;
  cin >> mc;
  initial_T = 1.0;
  final_T = 3.0;
  T_step = 0.1;
  int **spin_matrix = (int**) matrix(n, n, sizeof(int));
  long startpoint = -1;
  clock_t start, finish;
  double time_spent;
  start = clock();
  for(double T= initial_T; T <= final_T; T+=T_step){
    E = M = 0;
    for(int de=-8; de <= 8; de++) w[de+8] = 0;
    for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/T);
    for (int i = 0; i < 5; i++) average[i] = 0.;
    initialize(n,  T, spin_matrix, E, M);
    for(int cycles=1; cycles <= mc; cycles++){
      metropolis(n, startpoint, spin_matrix, E, M, w);
      average[0] += E;
      average[1] += E*E;
      average[2] += M;
      average[3] += M*M;
      average[4] += fabs(M);
  //    cout << "loop" << cycles << endl; //checing if the loop works
     }
     writingfunc(n, mc, T, average);
  }
  finish = clock();
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
  cout << "Time spent: " << time_spent << endl; //time spent on the loops and writing to file
  free_matrix((void**)spin_matrix); // freeing memory
  ofile.close();
  return 0;
}