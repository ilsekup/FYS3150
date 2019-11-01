#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "lib.h"
using namespace std;
ofstream ofile;

inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

void initialize(int n, double T, int **spin_matrix, double& E, double& M){
  for(int y=0; y < n; y++){ //columns
    for(int x=0; x < n; x++){ //rows
      if(T < 1.5){
        spin_matrix[x][y] = 1; //element is one if the temp is low
        }
      M+= (double)spin_matrix[x][y]; //updating Magnetization

    }
  }
  for (int y=0; y < n; y++){
    for(int x=0; x < n; x++){
      E -= (double)(spin_matrix[periodic(y, n, -1)][x] + spin_matrix[y][periodic(x, n, -1)]);

     }
  }
}

void metropolis(int n, long& startpoint, int **spin_matrix, double& E, double& M, double *w){
  for(int y = 0; y < n; y++){
    for(int x=0; x < n; x++){
      int ix = (int)(ran1(&startpoint)*(double)n);
      int iy = (int)(ran1(&startpoint)*(double)n);
      int dE = 2*spin_matrix[iy][ix]*(spin_matrix[iy][periodic(ix, n, -1)] + spin_matrix[periodic(iy, n,-1)][ix]
    + spin_matrix[iy][periodic(ix, n, 1)] + spin_matrix[periodic(iy, n, 1)][ix]);
    if(ran1(&startpoint) <= w[dE+8]){
      spin_matrix[iy][ix] *= -1;
      M += (double) 2*spin_matrix[iy][ix];
      E += (double) dE;
      }
    }
  }
}

void writingfunc(int n, int mc, double T, double *average)
{
double norm = 1/((double) (mc)); // divided by total number of cycles
double Eaverage = average[0]*norm;
double E2average = average[1]*norm;
double Maverage = average[2]*norm;
double M2average = average[3]*norm;
double Mabsaverage = average[4]*norm;
double Evariance = (E2average- Eaverage*Eaverage)/(n*n);
double Mvariance = (M2average - Maverage*Maverage)/(n*n);
double M2variance = (M2average - Mabsaverage*Mabsaverage)/(n*n);
double Mvarianceabs = (M2average - Mabsaverage*Mabsaverage)/(n*n);
ofile << setiosflags(ios::showpoint | ios::uppercase);
ofile << setw(15) << setprecision(8) << T;
ofile << setw(15) << setprecision(8) << Eaverage/(n*n);
ofile << setw(15) << setprecision(8) << Evariance/(T*T); // heat capacity
ofile << setw(15) << setprecision(8) << M2variance/T; //susceptibility
ofile << setw(15) << setprecision(8) << Mabsaverage/(n*n) << endl;
}

int main(int argc, char* argv[]){
  char *outfilename;
  long startpoint;
  int **spin_matrix, n, mc;
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
  spin_matrix = (int**) matrix(n, n, sizeof(int));
  startpoint = -1;
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

