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


int main(int argc, char* argv[])
{
  //initial setup
    char *outfilename;
    int n,t_steps;
    high_resolution_clock::time_point t1, t2;
    duration<double, ratio<1,1>> time_spent;

    omp_set_num_threads(4);
    int thread_num = omp_get_max_threads();
    cout << "  The number of processors available = " << omp_get_num_procs() << endl ;
    cout << "  The number of threads available    = " << thread_num <<  endl;

    if (argc<3){
      cout << "Use command line arguemnts: n (dimension of matrix), t_stop" << endl;
      return 0;
      }

    n = atof(argv[1]);
    t_steps = atof(argv[2]);
    double dx = 1/(double) (n+1);
    double dt = 0.4*dx*dx;

    t1 = high_resolution_clock::now();
    mat u = implicit2DMP(n,dt,t_steps,return_zero);
    t2 = high_resolution_clock::now();
    time_spent = t2-t1;
    cout << " Time_spent implicit 2D with OpenMP = " << \
    time_spent.count() << " seconds for n = " << n << " and t_steps = " << t_steps << endl;

    t1 = high_resolution_clock::now();
    mat u2 = implicit2D(n,dt,t_steps,return_zero);
    t2 = high_resolution_clock::now();
    time_spent = t2-t1;
    cout << " Time_spent implicit 2D = " << \
    time_spent.count() << " seconds for n = " << n << " and t_steps = " << t_steps << endl;

    return 0;
}
