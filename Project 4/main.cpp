#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "isingmodel.h"
#include "lib.h"

//must be compiled with lib.cpp and isingmodel.cpp and must be run with name of an outputfile

using namespace std;
int main(int argc, char* argv[]){
  char *outfilename;
  int n, Tlone, counter;
  double mc;
  double w[17], average[5], initial_T, final_T, E, M, T_step, energy[3];
  outfilename = argv[1];
  ofstream ofile2;
  ofstream ofile3;
  ofile2.open("mcsweeps.txt");
  ofile3.open("histdata.txt");
  cout << "Insert L, dimension of matrix: " << endl;
  cin >> n;
  cout << "Number of monte carlo cycles or max # sweeps" << endl;
  cin >> mc;
  cout << "Insert constant temperature T = 1.0 or T = 2.4" << endl;
  cin >> Tlone;
  int **spin_matrix = (int**) matrix(n, n, sizeof(int));
  long startpoint = -1;
  clock_t start, finish;
  double time_spent;
  long seed = time(NULL);
  start = clock();

  bool timeplot = true; //this loop is only needed for exercise c, for the time dependence
  if(timeplot == true)
 {
  for(int sweep = 1; sweep < mc; sweep = sweep + 10)
    {
    double E2, M2;
    E2 = M2 = 0;
    counter = 0;
    for(int de2=-8; de2 <= 8; de2++) w[de2+8] = 0;
    for(int de2 =-8; de2 <= 8; de2+=4)w[de2+8] = exp(-de2/Tlone);
    for (int i = 0; i < 5; i++) energy[i] = 0.;
    initialize(n,  Tlone, spin_matrix, E2, M2, seed, true);  //false if we want random spins, true for ordered
    for(int cycles= 1; cycles <= sweep; cycles++)
        {
      metropolis(n, startpoint, spin_matrix, E2, M2, w, counter);
      energy[0] += E2;
      energy[1] += E2*E2;
      energy[3] += fabs(M2);
        }
    writingfunc2(n, sweep, Tlone, energy, counter, ofile2);
    }
 }

E = M = 0;
for(int de=-8; de <= 8; de++) w[de+8] = 0;
for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/Tlone);
for (int i = 0; i < 5; i++) average[i] = 0;
initialize(n,  Tlone, spin_matrix, E, M, seed, false);  //false if we want random spins, true for ordered
  for(int cycles=1; cycles <= mc; cycles++){
    metropolis(n, startpoint, spin_matrix, E, M, w, counter);
    double energy = E;
    average[0] += E;
    average[1] += E*E;
    writingfunc3(n, mc, Tlone, energy, average, cycles, ofile3);
    }

  finish = clock();
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
  cout << "Time spent: " << time_spent << endl; //time spent on the loops and writing to file

  free_matrix((void**)spin_matrix); // freeing memory
  ofile2.close();
  ofile3.close();

  return 0;
}
