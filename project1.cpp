#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>

std::ofstream ofile;
using namespace std;

double f(double x){
  return 100*exp(-10*x);
}

double u(double x){
  return 1 - (1- exp(-10))*x - exp(-10*x);
}

void LUdecomposition(double *b_tilde, double n){
  using namespace arma;
  mat A = zeros<mat>(n,n);
  vec b_tilde2(n);
  for (int i = 0; i < n - 1; i++){ //making the tridiagonal matrix
    A(i,i) = 2;
    A(i,i+1) = -1;
    A(i+1,i) =-1;
    b_tilde2[i] = b_tilde[i];
  }
  A(n-1,n-1) = 2; // setting the last corner element to -2;
  mat L, U;
  lu(L,U,A); // LU-decomposition function from armadillo
  vec z = solve(L,b_tilde2); // solve for z = U*v, Lz = b_tilde
  vec v = solve(U,z); //solve for v
  //v.print("v=");
}

void solver_Thompson(double *a, double *b, double *c, double *b_tilde, double *&v, int n){
    // Solves the matrix problem (for tridiagonal matrix)
    for(int i = 2; i < n-1; i++){
        b[i] -= c[i-1]*a[i-1]/b[i-1]; // 3 FLOPS
        b_tilde[i] -= b_tilde[i-1]*a[i-1]/b[i-1]; //3 FLOPS
    }
    v[n-2] = b_tilde[n-2]/b[n-2];
    for(int i = n-3; i>0; i--){
        v[i] = (b_tilde[i]-(v[i+1]/c[i]))/b[i];// 3 FLOPS
    }
    // sum 9N FLOP
}


void solver_Thompson_quick(double *b, double *b_tilde, double *&v, int n){
    // Solves the matrix problem in the simplified way
    // Assumes a = -1, b = 2 and c = -1
    // Requires that b is precalculated by b_i = (i+1)/i
    // Forwards substitution
    for(int i = 2; i < n-1; i++){
        b_tilde[i] += b_tilde[i-1]/b[i-1]; // 2 FLOPS
    }
    // Backwards substitution
    v[n-2] = b_tilde[n-2]/b[n-2];
    for(int i = n-3; i>0; i--){
        v[i] = (b_tilde[i]+v[i+1])/b[i]; // 2 FLOPS

    }
    // sum 4N FLOPS
}

void print_file(double error, double t, double *v, double *x, string outfilename, int k){
    // prints the results in a textfile
    // includes the numerical solution v if n <= 1e4
    ofile.open(outfilename);
    ofile << "Max error = " << error << endl;
    ofile << "Time used = " << t << endl;
    if (k<=1e4){
        ofile << "   x  ,   u   ,   v" << endl;
        ofile << "------------------------------" << endl;
        for(int i = 0; i < k+2; i++){
        ofile << showpoint << setprecision(6) << setw(6) << x[i] << " , " << u(x[i]) << " , " << v[i] <<endl;
    }
    }
    ofile.close();
}

double find_max_error(double *u, double *v, int n){
    // Returns log_10 of the maximum relative error
    // u is the exact solution and v is the calculated, n is the number of points in u and v
    // exludes u[0] and u[n-1] as these are 0
    double rel_error;
    double max_error = log10(fabs((v[1] - u[1])/u[1]));
    for(int i = 2; i<n-1; i++){
        rel_error = log10(fabs((v[i] - u[i])/u[i])); //logarithm of relative error
        if (rel_error>max_error){
        max_error = rel_error;
        }
    }
    return max_error;


}

int main(int argc, char *argv[]){
  //initialize vectors for a, b, c, b_tilde, x and v
  int n = atof(argv[1]);//number of points we calculate on
  int N = n+2;//actual number of elements (includes x_0 and x_n+1)
  clock_t start, finish;
  double *v = new double[N];
  double *b = new double[N];
  double *b_tilde = new double[N];
  double *c = new double[N];
  double *a = new double[N];
  double *x = new double[N];
  double *exact = new double[N];

  double h = 1.0/(n+1);

  //initialising first elements in the vectors
  x[0] = 0, x[n+1] = 1;
  exact[0] = exact[n+1] = 0;
  double hh = h*h; //do this outside loop to minimize flops in loop

  //doing a loop for initializing
  std::cout << "Initializing vectors, and calculating exact vectors" << std::endl;
  for(int i = 1; i < n+1; i++){
    a[i] = -1, c[i] =-1, b[i] = 2; //initialising vecotors a, b, c
    x[i] = i*h; // making vector x
    b_tilde[i] = hh*f(x[i]); // making vector b_tilde
    exact[i] = u(x[i]); // finding exact solution
  }

  start = clock();
  // calculate the solution with the 'slow' algorithm
  solver_Thompson(a,b,c,b_tilde,v,N);
  finish = clock();
  // free up space since I don't need a or c any more
  delete [] a;
  delete [] c;

  // find the error in our numerical solution
  std::cout<<"Calculating the error between calculated and expected answer \n";
  double max_error = find_max_error(exact,v,N);

  // calculate the time it took
  double time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );

  cout << "Max relative error = " << max_error << endl;
  cout <<"Time elapsed long version =  " << time_spent <<" seconds" <<endl;

  // save results to file
  string outfilename;
  string number = argv[1];
  outfilename = "output_"+number+".data";

  print_file(max_error,time_spent,v,x, outfilename,n);

  // restarts, with the faster algorithm
  cout << "Restart, test with quicker algorithm \n\n"<<endl;
  for(int i = 1; i <= n; i++){
    b[i] = (i+1)/(double)i; //precalculated values for b
    b_tilde[i] = hh*f(x[i]); // resets b_tilde
  }
  start = clock();
  // calculate the solution with the 'slow' algorithm
  solver_Thompson_quick(b,b_tilde,v,N);
  finish = clock();
  // find the error in our numerical solution, this should be (about) the same
  std::cout<<"Calculating the error between calculated and expected answer \n";
  max_error = find_max_error(exact,v,N);

  // calculate the time it took
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );

  cout << "Max relative error = " << max_error << endl;
  cout <<"Time elapsed long version =  " << time_spent <<" seconds" <<endl;

  // save results to file
  outfilename = "output_quick_"+number+".data";

  print_file(max_error,time_spent,v,x,outfilename,n);

  cout<<"Calculating v using the LU-decomposition\n";
  //timing the LU decomposition
  start = clock();
  LUdecomposition(b_tilde, N);
  finish = clock();
  time_spent = ((double)(finish - start)/ CLOCKS_PER_SEC);
  cout <<"Time elapsed LU = " << time_spent << " seconds" <<endl;
  return 0;
}
