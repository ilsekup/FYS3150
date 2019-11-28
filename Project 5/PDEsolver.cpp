#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "lib.h"
#include <armadillo>
#include "PDEsolver.h"

using namespace std;
using namespace arma;

ofstream ofile("explicit.txt", ios::out);
ofstream ofile2("implicit.txt", ios::out);
ofstream ofile3("explicit2D.txt", ios::out);




void writingfunc1D(int n, vec u_t, ostream& ofile)
{
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for (int k = 0; k <= n; k++)
    {
      ofile << setw(15) << setprecision(8) << u_t(k); // writes u for all x
    }
  ofile << endl;
}



void writingfunc2D(int n, mat u_t, ostream& ofile3)
{
  ofile3 << setiosflags(ios::showpoint | ios::uppercase);

  for (int i = 0; i <= n; i++)
    {
      ofile3 << endl;
      for (int k = 0; k <= n; k++)
      {
      ofile3 << setw(15) << setprecision(8) << u_t(i,k); // writes u for all x
      }
    }
  ofile3 << endl;
}




//explicit 1D
void explicitsch1D(int n, int t_steps)
{
  // setting up steplength in x and t
  double dx = 1.0 / (n+1);
  double dt = 1.0 /(t_steps);
  double alpha = dt / ( dx * dx);

  // setting up vectors and initial/bouandry conditions
  vec u_xx = zeros<vec>(n+1);
  vec u_t = zeros<vec>(n+1);
  u_xx(0) = u_t(0) = 0.0;
  u_xx(n) = u_t(n) = 1.0;

  writingfunc1D(n, u_t, ofile); // to write first line at i,j = 0
  for (int j = 1; j < t_steps + 1; j++) // iterating over temperatures
    {
      for (int i = 1; i < n; i++) // iterating over x-position
      {
      u_t(i) = alpha*( u_xx(i-1) + u_xx(i+1) ) + (1 - 2*alpha) * u_xx(i);
      u_xx(i) = u_t(i);
      }
      writingfunc1D(n, u_t, ofile);
  }
  ofile.close();
}

// explicit 2D
void explicitsch2D(int n, int t_steps)
{
  double dt = 1.0 /(t_steps);
  // setting up steplength in x = y and t
  double dx = 1.0 / (n+1); //dy = dx , alpha is the same for both as well
  double alpha = dt / ( dx * dx);

  // setting up vectors/matrices and initial/bouandry conditions
  mat u_yx = zeros<mat>(n+1,n+1);
  mat u_t = zeros<mat>(n+1,n+1);


  //these conditions imply that the left and right side of the lattice is 1, 0 everywhere else
  for(int i = 0; i < n + 1 ; i++)
  {
    u_yx(i,0) = 1.0;
    u_yx(i,n) = 1.0;
    u_yx(0,i) = 0.0;
    u_yx(n,i) = 0.0;
  }
  u_yx(n,0) = 1;
  u_yx(0,n) = 1;
  u_yx(n,n) = 1;
  u_yx(0,0) = 1;

u_t = u_yx; //setting matrices equal to each other so they have the same initial/bouandry conditions

ofile3 << setiosflags(ios::showpoint | ios::uppercase);
ofile3 << n+1 << endl;

writingfunc2D(n, u_t, ofile3); // to write first line at i,j = 0

for (int j = 1; j < t_steps; j++) // iterating over temperatures
  {
    for (int i = 1; i < n; i++) // iterating over x-position
    {
        for(int k = 1; k < n; k++) //iteratinv over y-position
        {
        u_t(i,k) = alpha*(u_yx(i-1,k) + u_yx(i+1,k) + u_yx(i,k+1) + u_yx(i,k-1)) + (1 - 4*alpha) * u_yx(i,k);
        u_yx(i,k) = u_t(i,k);
        }
    }

    writingfunc2D(n, u_t, ofile3);
  }
  ofile3.close();
}

void writingfunc2(int n, double *u, ostream& ofile){
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for (int k = 0; k < n; k++)
  {
    ofile << setw(15) << setprecision(8) << u[k]; // writes u for all x
  }
  ofile << endl;
}

void solver_Thomas(double *a, double *b, double *c, double *vold, double *vnew, int n){
    // Solves the matrix problem (for tridiagonal matrix)
    double *bt = new double[n+1];
    double *voldt = new double[n+1];
    bt[1] = b[1];
    voldt[1] = vold[1];
    for(int i = 2; i < n+1; i++){
        bt[i] = b[i] - c[i-1]*a[i-1]/bt[i-1]; // 3 FLOPS
        voldt[i] = vold[i]- voldt[i-1]*a[i-1]/bt[i-1]; // 3 FLOPS
    }
    vnew[n] = voldt[n]/bt[n]; // 1 FLOP
    for(int i = n-1; i>0; i--){
      vnew[i] = (voldt[i]-vnew[i+1]*c[i])/bt[i];// 3 FLOPS
    }
    // sum 9N + 1 FLOP
}

void solver_Thomas2(double* a, double* b, double* c, double* v, double* u, int n){
    // Solves the matrix problem Av=u (for tridiagonal matrix A with a, b and c along the diagonal)
    // Vectors b,v and u must be n long, while a and c are n-1 long (as they are the off diagonal elements)
    // Overwrites b, u and v

    for(int i = 1; i < n; i++){
      b[i] -= a[i-1]*c[i-1]/b[i-1];
      u[i] -= u[i-1]*a[i-1]/b[i-1];
    }

    v[n-1] = u[n-1]/b[n-1];
    for(int i = n-2; i>=0; i--){
      v[i] = (u[i]-v[i+1]*c[i])/b[i];
    }
}

double *implicit(int n, double t_steps){
  double dx = 1/double(n);
  double dt = 0.5*dx*dx;
  double alpha = dt/(dx*dx);
  double *vold = new double[n+2];
  double *vnew = new double[n+2];
  double *vinitial = new double[n+2];
  double *a = new double[n];
  double *b = new double[n];
  double *c = new double[n];
  //initial and boundary conditions
  for(int j = 0; j< n; j++){
    vold[j] = 0;
  }
  for(int i= 0; i <= n; i++){
    vnew[0] = vold[0] = 0;
    vnew[n] = vold[n] = 1;
    c[i] = a[i] = -alpha;
    b[i] = 1 + alpha*2;
    vinitial[i] = vold[i];
  }
  for(int t = 1; t < t_steps; t++){
    writingfunc2(n, vnew, ofile2);
    solver_Thomas(a, b, c, vold, vnew, n);
    vnew[0] = 0;
    vnew[n] = 1;
    for(int i = 0; i <= n; i++){
      vold[i] = vnew[i];

     }
  }
  ofile.close();
  return vold;
}

// One-dimensional Crank–Nicolson
void CN(int n, double dt, int t_steps, ostream& ofile){
  double dx = 1/(double) (n+1);
  double alpha = dt / (double) (dx*dx);
  double *v = new double[n+2];
  double *v_tilde = new double[n+2];
  double *a = new double[n+1];
  double *b = new double[n+2];
  double *c = new double[n+1];
  double b_start = 2 + 2*alpha;;
  double _2minus2alpha = 2.0 - 2.0*alpha; // precalculated for speed

  // Initialize values
  for(int i=1; i<n+1; i++){
    b[i] = b_start;
    v[i] = 0;
  }
  b[0] = b[n+1] = 1;
  v[0] = 0;
  v[n+1] = 1;

  for(int i=1; i<n; i++) a[i] = c[i] = -alpha;
  a[n] = c[0] = 0;
  a[0] = c[n] = -alpha;

  writingfunc2(n+2, v, ofile2);

  // Solve for time
  for(int t=0; t<t_steps; t++){
    // Find v_tilde
    v_tilde[0] = v[0];
    v_tilde[n+1] = v[n+1];
    for(int j=1; j<=n; j++){
      v_tilde[j] = alpha*v[j-1] + _2minus2alpha * v[j] + alpha*v[j+1];
      //Reset b
      b[j] = b_start;
    }
    b[0] = b[n+1] = 1;

    solver_Thomas2(a,b,c,v,v_tilde,n+2);
    writingfunc2(n+2, v, ofile2);
  }
}
