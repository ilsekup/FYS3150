#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "PDEsolver.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  //setup for writing in to file
  char *outfilename;
  outfilename = argv[1];

  //choosing n ad t steps
  int n, t; // number of steps in x and t respectively.
  double t_stop;
  cout << "n = ";
  cin >> n;
  cout << "t_stop = ";
  cin >> t_stop;

  double dx = 1/(double) (n+1);
  double dt = 0.4*dx*dx;
  t = t_stop/dt;

  ofstream ofile(outfilename, ios::out);
  //calling function which also does the writing into file
  explicitsch1D(n,dt,t);
  implicit(n, dt, t);
  explicitsch2D(n,dt,t);
  CN(n,dt,t,ofile);
  int it = implicit2D(n, dt, t);
  cout << it << endl;
  ofile.close();

  ofstream ofile_info("runinfo.txt", ios::out);
  ofile_info << "dt = " << dt << " n = " << n  << " t_stop = " << t_stop << endl;
  ofile_info.close();
  return 0;
}
