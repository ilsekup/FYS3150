#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "PDEsolver.h"
#include <armadillo>

using namespace std;
using namespace arma;

// Filler function to solve laplace with a poisson solver
double return_zero(double t, double x, double y){
  return 0.0;
}

mat set_initial_conditions(int n){
  mat u = zeros<mat>(n+2, n+2); //timesteps before

  // Setting initial and boundary conditions (most is just 0)
  for(int j=0; j<n+2;j++){
    u(0, j) = 1.0;
    u(n+1, j) = 1.0;
  }
  return u;
}
int main(int argc, char* argv[])
{
  //setup for writing in to file
  char *outfilename;
  outfilename = argv[1];

  //choosing n and t steps
  int n, t; // number of steps in x and t respectively.
  double t_stop;
  cout << "n = ";
  cin >> n;
  cout << "t_stop = ";
  cin >> t_stop;

  int nx,ny;
  nx = ny = n;
  double dx = 1/(double) (n+1);
  double dt = 0.4*dx*dx;
  t = t_stop/dt;

  //calling function which also does the writing into file
  explicitsch1D(n,dt,t);
  implicit(n, dt, t);
  mat u_t = explicitsch2D(n,dt,t);

  ofstream ofile(outfilename, ios::out);
  CN(n,dt,t,ofile);
  ofile.close();

  ofstream ofile2("implicit2D.txt", ios::out);
  mat u0 = set_initial_conditions(n);
  mat u = implicit2D(nx, ny, dt, t, u0, return_zero, ofile2);
  ofile2.close();

  ofstream ofile_info("runinfo.txt", ios::out);
  ofile_info << "dt = " << dt << " nx = " << n << " ny = " << n  << " t_stop = " << t_stop << endl;
  ofile_info.close();
  return 0;
}
