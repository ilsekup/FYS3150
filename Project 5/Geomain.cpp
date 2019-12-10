#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "PDEsolver.h"
#include <armadillo>

using namespace std;
using namespace arma;

double no_heat(double t, double x, double y){
  return 0;
}

// Returns the heat production as a function of position and time
// (Neglects time evolution)
double heat(double t, double x, double y){

  double W1 = 8064.0/1573.0;
  double W2 = 2016.0/1573.0;
  double W3 = 288.0/1573.0;
  double val=0.0;

  if (x<1.0/6.0){
    val=W1;
  } else if (x<1.0/3.0){
    val=W2;
  } else {
    val=W3;
  }
  return val;
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

  double dx = 1/(double) (n+1);
  double dt = 0.4*dx*dx;
  t = t_stop/dt;

  ofstream ofile(outfilename, ios::out);
  //calling function which also does the writing into file
  mat u = implicit2D(n, dt, t, heat);
  ofile.close();

  ofstream ofile_info("runinfo.txt", ios::out);
  ofile_info << "dt = " << dt << " n = " << n  << " t_stop = " << t_stop << endl;
  ofile_info.close();
  return 0;
}
