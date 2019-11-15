#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "isingmodel.h"
#include "lib.h"

//must be compiled with lib.cpp and must be run with name of an outputfile

using namespace std;
int main(int argc, char* argv[]){
  char *outfilename;
  int n, Tlone;
  double mc;
  double w[17], average[5], initial_T, final_T, E, M, T_step, energy[3];
;
  outfilename = argv[1];
  ofstream ofile;
  ofile.open(outfilename);
  ofile2.open("temp.txt");
  ofile3.open("histdata.txt");
  cout << "Insert L, dimension of matrix: " << endl;
  cin >> n;
  cout << "Number of monte carlo cycles" << endl;
  cin >> mc;
  cout << "Insert constant temperature T = 1.0 or T = 2.4" << endl;
  cin >> Tlone;
  initial_T = 1.0;
  final_T = 3.0;
  T_step = 0.1;
  int **spin_matrix = (int**) matrix(n, n, sizeof(int));
  long startpoint = -1;
  clock_t start, finish;
  double time_spent;
  long seed = time(NULL);
  start = clock();
  for(double T= initial_T; T <= final_T; T+=T_step){
    E = M = 0;
    for(int de=-8; de <= 8; de++) w[de+8] = 0;
    for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/T);
    for (int i = 0; i < 5; i++) average[i] = 0.;
    initialize(n,  T, spin_matrix, E, M, seed, false);  //false if we want random spins, true for ordered
    for(int cycles=1; cycles <= mc; cycles++){
      metropolis(n, startpoint, spin_matrix, E, M, w);
      average[0] += E;
      average[1] += E*E;
      average[2] += M;
      average[3] += M*M;
      average[4] += fabs(M);
     }
     writingfunc(n, mc, T, average,ofile);
  }
  bool timeplot = false; //this loop is only needed for exercise c, for the time dependence
  if(timeplot == true){
  for(int sweep = 1; sweep < 10000; sweep = sweep + 10){
    double E2, M2;
    E2 = M2 = 0;
    for(int de2=-8; de2 <= 8; de2++) w[de2+8] = 0;
    for(int de2 =-8; de2 <= 8; de2+=4)w[de2+8] = exp(-de2/Tlone);
    for (int i = 0; i < 5; i++) energy[i] = 0.;
    initialize(n,  Tlone, spin_matrix, E2, M2, seed, true);  //false if we want random spins, true for ordered
    for(int cycles= 1; cycles <= sweep; cycles++){
      metropolis(n, startpoint, spin_matrix, E2, M2, w);
      energy[0] += E2;
      energy[1] += E2*E2;
      energy[3] += fabs(M2);

    }
    writingfunc2(n, sweep, Tlone, energy);
  }
}
E = M = 0;
for(int de=-8; de <= 8; de++) w[de+8] = 0;
for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/Tlone);
for (int i = 0; i < 5; i++) average[i] = 0;
initialize(n,  Tlone, spin_matrix, E, M, seed, false);  //false if we want random spins, true for ordered
  for(int cycles=1; cycles <= mc; cycles++){
    metropolis(n, startpoint, spin_matrix, E, M, w);
    double energy = E;
    average[0] += E;
    average[1] += E*E;
    writingfunc3(n, mc, Tlone, energy, average, cycles);
    }
  finish = clock();
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
  cout << "Time spent: " << time_spent << endl; //time spent on the loops and writing to file
  free_matrix((void**)spin_matrix); // freeing memory
  ofile2.close();
  ofile.close();
  ofile3.close();

  return 0;
}