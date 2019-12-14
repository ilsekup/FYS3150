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

void writingfunc2D(int nx, int ny, mat u_t, ostream& ofile3)
{
  ofile3 << setiosflags(ios::showpoint | ios::uppercase);

  for (int i = 0; i < nx; i++)
    {
      ofile3 << endl;
      for (int k = 0; k < ny; k++)
      {
      ofile3 << setw(15) << setprecision(8) << u_t(i,k); // writes u for all x
      }
    }
  ofile3 << endl;
}




//explicit 1D
void explicitsch1D(int n, double dt, int t_steps)
{
  // setting up steplength in x and t
  double dx = 1.0 / (n+1);
  double alpha = dt / ( dx * dx);
  ofstream ofile("explicit.txt", ios::out);

  // setting up vectors and initial/bouandry conditions
  vec u_xx = zeros<vec>(n+2);
  vec u_t = zeros<vec>(n+2);
  u_xx(0) = u_t(0) = 0.0;
  u_xx(n+1) = u_t(n+1) = 1.0;

  u_t.t().raw_print(ofile); // to write first line at i,j = 0
  for (int j = 1; j < t_steps + 1; j++) // iterating over temperatures
    {
      for (int i = 1; i < n+1; i++) // iterating over x-position
      {
      u_t(i) = alpha*( u_xx(i-1) + u_xx(i+1) ) + (1 - 2*alpha) * u_xx(i);
      u_xx(i) = u_t(i);
      }
      u_t.t().raw_print(ofile);
  }
  ofile.close();
}

// explicit 2D
mat explicitsch2D(int n, double dt, int t_steps)
{
  // setting up steplength in x = y and t
  double dx = 1.0 / (n+1); //dy = dx , alpha is the same for both as well
  double alpha = dt / ( dx * dx);
  ofstream ofile("explicit2d.txt", ios::out);

  // setting up vectors/matrices and initial/bouandry conditions
  mat u_yx = zeros<mat>(n+2,n+2);
  mat u_t = zeros<mat>(n+2,n+2);


  //these conditions imply that the left and right side of the lattice is 1, 0 everywhere else
  for(int i = 0; i < n + 1 ; i++)
  {
    u_yx(i,0) = 1.0;
    u_yx(i,n+1) = 1.0;
    u_yx(0,i) = 0.0;
    u_yx(n+1,i) = 0.0;
  }
  u_yx(n+1,0) = 1;
  u_yx(0,n+1) = 1;
  u_yx(n+1,n+1) = 1;
  u_yx(0,0) = 1;

  u_t = u_yx; //setting matrices equal to each other so they have the same initial/bouandry conditions

  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << n+2 << endl;

  writingfunc2D(n+2,n+2, u_t, ofile); // to write first line at i,j = 0

  for (int j = 1; j < t_steps; j++) // iterating over time
  {
    for (int i = 1; i < n+1; i++) // iterating over x-position
    {
        for(int k = 1; k < n+1; k++) //iteratinv over y-position
        {
        u_t(i,k) = alpha*(u_yx(i-1,k) + u_yx(i+1,k) + u_yx(i,k+1) + u_yx(i,k-1)) + (1 - 4*alpha) * u_yx(i,k);
        u_yx(i,k) = u_t(i,k);
        }
    }

    writingfunc2D(n+2,n+2, u_t, ofile);
    //u_t.print(" ");
  }
  ofile.close();
  //u_yx.print(" yx");
  return u_yx;
}

void writingfunc(int n, double *u, ostream& ofile){
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for (int k = 0; k < n; k++)
  {
    ofile << setw(15) << setprecision(8) << u[k]; // writes u for all x
  }
  ofile << endl;
}

void solver_Thomas(double* a, double* b, double* c, double* u, double* v, int n){
    // Solves the matrix problem Av=u (for tridiagonal matrix A with a, b and c along the diagonal)
    // Vectors b,v and u must be n long, while a and c are n-1 long (as they are the off diagonal elements)
    // Overwrites b, u and v

    // Forwards substitution
    for(int i = 1; i < n; i++){
      b[i] -= a[i-1]*c[i-1]/b[i-1];
      u[i] -= u[i-1]*a[i-1]/b[i-1];
    }
    // Backwards substitution
    v[n-1] = u[n-1]/b[n-1];
    for(int i = n-2; i>=0; i--){
      v[i] = (u[i]-v[i+1]*c[i])/b[i];
    }
}

double *implicit(int n, double dt, double t_steps){
  double dx = 1/(double)(n+1);
  double alpha = dt/(dx*dx);
  double *vold = new double[n+2];
  double *vnew = new double[n+2];
  ofstream ofile2("implicit.txt", ios::out);
  double *a = new double[n+1];
  double *b = new double[n+2];
  double *c = new double[n+1];
  //initial and boundary conditions
  for(int j = 0; j <= n; j++){
    vold[j] = vnew[j] = 0;
  }

  vnew[0] = vold[0] = 0;
  vnew[n+1] = vold[n+1] = 1;

  for(int i=1; i<=n; i++){
    c[i] = a[i] = -alpha;
    b[i] = 1 + alpha*2;
  }
  a[n] = c[0] = 0;
  a[0] = c[n] = -alpha;
  b[0] = b[n+1] = 1;

  for(int t = 1; t < t_steps; t++){
    writingfunc(n+2, vnew, ofile2);
    solver_Thomas(a, b, c, vold, vnew, n+2);
    vnew[0] = 0;
    vnew[n+1] = 1;
    for(int i = 0; i <= n; i++){
      vold[i] = vnew[i];
      b[i] = 1 + alpha*2;
    }
    b[0] = b[n+1] = 1;
  }
  ofile2.close();
  return vold;
}

// One-dimensional Crank–Nicolson
void CN(int n, double dt, int t_steps, ostream& ofile){
  // n is the number of point exluding boundaries
  // We therefor have a total of n+2 points
  // results are printed to ofile
  double dx = 1/(double) (n+1);
  double alpha = dt / (double) (dx*dx);
  // Allocate required vectors
  double *v = new double[n+2];
  double *v_tilde = new double[n+2];
  double *a = new double[n+1]; // vecor a (and c) are one shorter, as they are off-diagonal
  double *b = new double[n+2];
  double *c = new double[n+1];
  double b_start = 2 + 2*alpha;;
  double _2minus2alpha = 2.0 - 2.0*alpha; // precalculated for speed

  for(int i=1; i<n+1; i++){
    // Set initial conditions and b vector
    b[i] = b_start;
    v[i] = 0;
  }
  // Set up boundary conditions
  b[0] = b[n+1] = 1;
  v[0] = 0;
  v[n+1] = 1;

  // Set up a and c vector
  for(int i=1; i<n; i++) a[i] = c[i] = -alpha;
  a[n] = c[0] = 0;
  a[0] = c[n] = -alpha;

  writingfunc(n+2, v, ofile);

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

    // Second step of Cranck-Nicolson
    solver_Thomas(a,b,c,v_tilde,v,n+2);
    // Save reusults
    writingfunc(n+2, v, ofile);
  }
}

// Solves 2D heat equation (with extra term) with implicit scheme and Jacobi solver
// nx and ny is the number of points in x and y directions respectively
// dt is the timestep, and t_step the total number of steps
// u is the full initial conditons this will in practice also give boundary conditions, as these are not changed
// (This also means that boundary conditions must be fixed values for this function to work)
// f is the functions of t,x,y (in that order) that gives the extra term in the heatequation (heatproduction, set to 0 for standard)
// ofils is the ostream where results are printed
mat implicit2D(int nx, int ny, double dt, int t_steps, mat &u, double (*f)(double,double,double),ostream& ofile){
  double h = 1/(double) (nx+1);
  double tol = 1e-8;
  double maxiterations = 10000;
  mat rho_tilde = zeros<mat>(nx+2,ny+2);
  writingfunc2D(nx+2, ny+2, u, ofile); // to write first line at i,j = 0
  for(int t= 0; t < t_steps; t++){
    for(int i=0;i<nx+2;i++){
      for(int j=0;j<ny+2;j++){
        rho_tilde(i,j) = f(t*dt,i*h,j*h);
      }
    }
    u = jacobi(nx, ny, dt, h, u,rho_tilde,tol);
    writingfunc2D(nx+2, ny+2, u, ofile); // to write first line at i,j = 0
  }
  return u;
}

// The jacobi solver for the implicit 2d scheme
mat jacobi(int nx, int ny, double dt, double h, mat &u, mat &rho_tilde, double tol){
    double alpha = dt/double(h*h);
    mat uold = zeros<mat>(nx+2, ny+2); //timesteps before
    mat uprev = zeros<mat>(nx+2, ny+2); //timesteps before

    for(int i = 1; i < nx+1; i++){ //setting up an initial guess for the old timestep
        for(int j = 1; j < ny+1; j++){
          uold(i,j) = 1.0;
          uprev(i,j) = u(i,j);
          }
        }

    for(int k = 0; k < 10000; k++){
    for(int i=1; i < nx+1; i++){
      for(int j=1; j < ny+1; j++){
        u(i,j) = (alpha*(uold(i+1, j) + uold(i-1, j) + uold(i, j+1) + uold(i, j-1)) + uprev(i, j)+ rho_tilde(i,j)*dt)/(double)(1+4*alpha);
      }
    }
    double sum = 0.0;
    for(int i = 0; i < nx+2; i++){
      for(int j =0; j < ny+2; j++){
        sum += (uold(i,j)- u(i,j))*(uold(i,j)-u(i,j));
        uold(i,j) = u(i,j);
      }
    }
    if(sqrt(sum) < tol){
      return u;
    }
    }
    return u;
}