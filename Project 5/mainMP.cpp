#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "PDEsolver.h"
#include <armadillo>
#include "omp.h"
#include "chrono"

using namespace std;
using namespace arma;
using namespace std::chrono;



double return_zero(double t, double x, double y){
  return 0;
}


mat implicit2DMP(int n, double dt, int t_steps, double (*f)(double,double,double)){ //implicit with jacobi solver + OMP
  double dx = 1/(double) (n+1);
  double tol = 1e-8;
  double maxiterations = 10000;
  mat A = zeros<mat>(n+2,n+2);
  mat u = jacobi(n, dt, dx, A,f,tol, t_steps);
  return u;
}

mat jacobiMP(int n, double dt, double h, mat &u, double (*f)(double,double,double), double tol, int t_steps){ //+OMP
    double dx = 1/(double)(n+1);
    double alpha = dt/double(dx*dx);
    mat uold = zeros<mat>(n+2, n+2); //timesteps before
    mat rho_tilde = zeros<mat>(n+2,n+2);
    ofstream ofile("implicit2dMP.txt", ios::out);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << n+2 << endl;
    // Setting boundary conditions
    for(int i=0; i < n+2; i++){
      u(0, i) = 1.0;
      u(n+1, i) = 0.0;
      u(i, 0) = 1-i*dx;
      u(i, n+1) = 1-i*dx;
    }
    u(0, n+1) = 1.0;
    u(n+1, 0) = 0.0;
    for(int i = 0; i < n+2; i++){ //setting up an initial guess for the old timestep
        for(int j = 0; j < n+2; j++){
          uold(i,j) = 1.0;
          }
        }
    writingfunc2D(n+2, u, ofile); // to write first line at i,j = 0
    for(int t= 0; t < t_steps; t++){
      for(int i=0;i<n+2;i++){
        for(int j=0;j<n+2;j++){
          rho_tilde(i,j) = f(t*dt,i*h,j*h)*h*h/4;
        }
      }
      for(int k = 0; k < 10000; k++){
        int i,j;
        double sum = 0.0;
        #pragma omp parlell default(shared) private(i,j) reduction (+:sum)
        { //start of paralell region
          # pragma omp for
      for(int i=1; i < n+1; i++){
        for(int j=1; j < n+1; j++){
          u(i,j) = (uold(i+1, j) + uold(i-1, j) + uold(i, j+1) + u(i, j-1))/4.0 + rho_tilde(i,j);
        }
      }

      for(int i = 0; i < n+2; i++){
        for(int j =0; j < n+2; j++){
          sum += (uold(i,j)- u(i,j))*(uold(i,j)-u(i,j));
          uold(i,j) = u(i,j);
        }
      }
    }//end of parllel region
      writingfunc2D(n+2, u, ofile);
      //u.print("u= ");
      if(sqrt(sum) < tol){
        return u;
      }
    }// end of jacobi iterations

    }// end of time loop
    ofile.close();
    return u;
}//end of jacobi function

mat explicitsch2DMP(int n)
{
  // setting up steplength in x = y and t
  double dx = 1.0 / (n+1); //dy = dx , alpha is the same for both as well
  int t_steps = 3*n*n;  //assures that t >= 2 n^2
  double dt = 1.0/t_steps;
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

  writingfunc2D(n+2, u_t, ofile); // to write first line at i,j = 0
  int i,k;
  for (int j = 1; j < t_steps; j++) // iterating over time
  {
  #pragma omp parlell default(shared) private(i,k)
  { //start of paralell region
      # pragma omp for
    for (int i = 1; i < n+1; i++) // iterating over x-position
    {
        for(int k = 1; k < n+1; k++) //iteratinv over y-position
        {
        u_t(i,k) = alpha*(u_yx(i-1,k) + u_yx(i+1,k) + u_yx(i,k+1) + u_yx(i,k-1)) + (1 - 4*alpha) * u_yx(i,k);
        u_yx(i,k) = u_t(i,k);
        }
    }
  }
    writingfunc2D(n+2, u_t, ofile);
    //u_t.print(" ");
  }
  ofile.close();
  //u_yx.print(" yx");
  return u_yx;
}

int main(int argc, char* argv[])
{
  //initial setup
    char *outfilename;
    int n,t_steps;
    high_resolution_clock::time_point t1, t2;
    duration<double, ratio<1,1>> time_spent;

    int thread_num = omp_get_max_threads();
    cout << "  The number of processors available = " << omp_get_num_procs() << endl ;
    cout << "  The number of threads available    = " << thread_num <<  endl;

    ofstream ofile2("Time.txt", ios::out);
    ofile2 << setiosflags(ios::showpoint | ios::uppercase);
    for(n = 10; n <= 100; n += 5)
    {
    ofile2 << setw(15) << setprecision(8) << n;
    for (int p = 1; p <= 4; p++) // p = 1, 2, 3, 4
     {
      omp_set_num_threads(p);
      t1 = high_resolution_clock::now();
      mat u = explicitsch2DMP(n);
      t2 = high_resolution_clock::now();
      time_spent = t2-t1;
      cout << " Time spent on Explicit 2D with OpenMP = " << \
      time_spent.count() << " seconds for n = " << n  << " and " << p << " thread(s)  "<< endl;
      ofile2 << setw(15) << setprecision(8) << time_spent.count();
    }
  ofile2 << endl;
  }
  ofile2.close();
  return 0;
}
