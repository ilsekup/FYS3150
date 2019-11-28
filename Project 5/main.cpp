#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "PDEsolver.h"
#include "lib.h"
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
  cout << "n = ";
  cin >> n;
  cout << "t = ";
  cin >> t;

  double dx = 1/(double) (n+1);
  double dt = 0.5*dx*dx;
  ofstream ofile(outfilename, ios::out);
  //calling function which also does the writing into file
  explicitsch1D(n,t);
  implicit(n, t);
  explicitsch2D(n,t);
  CN(n,dt,t,ofile);
  return 0;
}
