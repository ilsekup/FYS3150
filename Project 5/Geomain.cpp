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
// (Neglects time evolution, for the first case)
double heat1(double t, double x, double y){
  double W1 = 8064.0/1573.0;
  double W2 = 2016.0/1573.0;
  double W3 = 288.0/1573.0;
  double val=0.0;

  if ((1-x)<1.0/6.0){
    val=W1;
  } else if ((1-x)<1.0/3.0){
    val=W2;
  } else {
    val=W3;
  }
  return val;
}

// Returns the heat production as a function of position and time
// (Neglects time evolution, for the second case)
double heat2(double t, double x, double y){
  double W1 = 8064.0/1573.0;
  double W2 = 2016.0/1573.0;
  double W3 = 288.0/1573.0;
  double W4 = 2880.0/1573.0;
  double val=0.0;

  if ((1-x)<1.0/6.0){
    val=W1;
  } else if ((1-x)<1.0/3.0){
    val=W2;
  } else {
    val=W3+W4;
  }
  return val;
}

// Returns the heat production as a function of position and time
// (with time evolution, for the third case)
double heat3(double t, double x, double y){
  double W1 = 8064.0/1573.0;
  double W2 = 2016.0/1573.0;
  double W3 = 288.0/1573.0;
  double W4 = 2880.0/1573.0;

  double halflife1 = 6.9931;
  double halflife2 = 21.902;
  double halflife3 = 1.9556; // confirmed

  double log2 = log(2.0);

  double tau1 = halflife1/log2;
  double tau2 = halflife3/log2;
  double tau3 = halflife3/log2;

  double val=0.0;

  if ((1-x)<1.0/6.0){
    val=W1;
  } else if ((1-x)<1.0/3.0){
    val=W2;
  } else {
    val= W3 + 0.4*W4*(exp(-t/tau1)+exp(-t/tau2)+0.5*exp(-t/tau3));
  }
  return val;
}

mat set_initial_conditions(int nx, int ny, double dx){
  mat u = zeros<mat>(nx+2, ny+2); //timesteps before

  // Setting initial and boundary conditions
  for(int i=0; i < nx+2; i++){
    for(int j=0; j<ny+2;j++){
      u(i, j) = 1-i*dx*0.8213;
    }
  }
  return u;
}

int main(int argc, char* argv[]){
  //choosing n and t steps
  int nx, ny, t; // number of steps in x and t respectively.
  double t_stop;
  cout << "n = ";
  cin >> nx;
  cout << "t_stop = ";
  cin >> t_stop;

  ny = round(1.25*nx);

  double dx = 1/(double) (nx+1);
  double dt = 0.4*dx*dx;
  t = t_stop/dt;
  mat u = set_initial_conditions(nx, ny, dx);

  cout << "Starting on case 1" << endl;
  ofstream ofile1("geosim1.txt", ios::out);
  mat u1 = implicit2D(nx, ny, dt, t, u, heat1, ofile1);
  ofile1.close();

  cout << "Starting on case 2" << endl;
  u = set_initial_conditions(nx, ny, dx);
  ofstream ofile2("geosim2.txt", ios::out);
  mat u2 = implicit2D(nx, ny, dt, t, u, heat2, ofile2);
  ofile2.close();

  // Simulate 3 Gyr with decay
  int t3;
  t3 = 3/dt;

  u = set_initial_conditions(nx, ny, dx);
  cout << "Starting on case 3 (3 Gyr simulation)" << endl;
  ofstream ofile3("geosim3.txt", ios::out);
  mat u3 = implicit2D(nx, ny, dt, t3, u, heat3, ofile3);
  ofile3.close();

  u = set_initial_conditions(nx, ny, dx);
  ofstream ofile_info("runinfo.txt", ios::out);
  ofile_info << "dt = " << dt << " nx = " << nx  << " ny = " << ny << " t_stop = " << t_stop << endl;
  ofile_info.close();
  return 0;
}
