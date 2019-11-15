#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "lib.h"
#include "isingmodel.h"
using namespace std;
ofstream ofile;
ofstream ofile2;


inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

void initialize(int n, double T, int **spin_matrix, double& E, double& M, long seed, bool ordered){
  srand(seed); //Randomize seed initialization
  for(int y=0; y < n; y++){ //columns
    for(int x=0; x < n; x++){ //rows
      if(ordered==true){
        spin_matrix[x][y] = 1;
      }
      else{
	       int randNum = rand() % 2; // Generate a random number between 0 and 1
         // cout << randNum << endl;
        if(randNum < 0.5){
          spin_matrix[x][y] = -1;
        }
        else{
          spin_matrix[x][y] = 1;
        }
      }
        // cout << spin_matrix[x][y] << endl;
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

void writingfunc2(int n, int mc, double T, double *energy)
{
  double norm = 1/((double) (mc)); // divided by total number of cycles
  double Eaverage = energy[0]*norm;
  double E2average = energy[1]*norm;
  double Mabsaverage = energy[3]*norm;
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  ofile2 << setw(15) << setprecision(8) << Eaverage/(n*n);
  ofile2 << setw(15) << setprecision(8) << Mabsaverage/(n*n) << endl;
  // for(int i; i < mc; i++){
  //   ofile2 << setw(15) << setprecision(8) << energy[i] / (n*n) << endl; // energy per particle
  // }

}
