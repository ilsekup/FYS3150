#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "isingmodel.h"
#include "isingmodel.cpp"

//must be compiled with lib.cpp and must be run with name of an outputfile
//writingfunction to make histogram
void writingfunc2(int n, int mc, double T, double E, double *average, int cycles)
{
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << cycles; // cycle number
  ofile << setw(15) << setprecision(8) << E / (n*n) << endl; // energy per particle
}

using namespace std;
int main(int argc, char* argv[])
{
  char *outfilename;
  int n;
  int mc;
  double w[17], T, E, M, average[5];
  outfilename = argv[1];
  ofile.open(outfilename);
  cout << "Insert L, dimension of matrix: " << endl;
  cin >> n;
  cout << "Insert constant temperature T = 1.0 or T = 2.4" << endl;
  cin >> T;
  cout << "Number of monte carlo cycles" << endl;
  cin >> mc;
  int **spin_matrix = (int**) matrix(n, n, sizeof(int));
  long startpoint = -1;
  E = M = 0;
  for(int de=-8; de <= 8; de++) w[de+8] = 0;
  for(int de =-8; de <= 8; de+=4)w[de+8] = exp(-de/T);
  for (int i = 0; i < 5; i++) average[i] = 0;
  initialize(n,  T, spin_matrix, E, M, false);  //false if we want random spins, true for ordered
    for(int cycles=1; cycles <= mc; cycles++){
      metropolis(n, startpoint, spin_matrix, E, M, w);
      double energy = E;
      average[0] += E;
      average[1] += E*E;
      writingfunc2(n, mc, T, energy, average, cycles);
    //    cout << "loop" << cycles << endl; //checing if the loop works
       }
    free_matrix((void**)spin_matrix); // freeing memory
    double norm = 1/((double) (mc)); // divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Evariance = (E2average - Eaverage*Eaverage)/(n*n);
    cout << "C_v = "<<Evariance / (T*T) << endl;
    cout << "< E > = " << Eaverage /(n*n) << endl;
    ofile.close();
    return 0;
}
