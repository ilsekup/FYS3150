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
using namespace arma;

double f(double x){
  return 100*exp(-10*x);
}

double u(double x){
  return 1 - (1- exp(-10))*x - exp(-10*x);
}

vec LUdecomposition(vec b_tilde, mat A){
  mat L, U;
  lu(L,U,A); // LU-decomposition function from armadillo
  vec z = solve(L,b_tilde); // solve for z = U*v, Lz = b_tilde
  vec v = solve(U,z); //solve for v
  return v;
}

void solver_Thomas(double *a, double *b, double *c, double *b_tilde, double *v, int n){
    // Solves the matrix problem (for tridiagonal matrix)
    // Overwrites b, b_tilde and v
    for(int i = 2; i < n-1; i++){
        b[i] -= c[i-1]*a[i-1]/b[i-1]; // 3 FLOPS
        b_tilde[i] -= b_tilde[i-1]*a[i-1]/b[i-1]; // 3 FLOPS
    }
    v[n-2] = b_tilde[n-2]/b[n-2]; // 1 FLOP
    for(int i = n-3; i>0; i--){
        v[i] = (b_tilde[i]-(v[i+1]/c[i]))/b[i];// 3 FLOPS
    }
    // sum 9N + 1 FLOP
}


void solver_Thomas_Specialized(double *b, double *b_tilde, double *v, int n){
    // Solves the matrix problem in the simplified way
    // Assumes a = -1, b = 2 and c = -1
    // Requires that b is precalculated by b_i = (i+1)/i
    // Forwards substitution
    // Overwrites b, b_tilde and v
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

void print_file(double error, double t, double *v, double *x, string outfilename, int k, string run_info_file){
    // prints the results in a textfile
    // prints the numerical solution to "outfilename" if n <= 1e4
    if (k<=1e4){

        ofile.open(outfilename);
        ofile << "Max error = " << error << endl;
        ofile << "Time used = " << t << endl;
        ofile << "   x  ,   u   ,   v" << endl;
        ofile << "------------------------------" << endl;
        for(int i = 0; i < k+2; i++){
          ofile << showpoint << setprecision(6) << setw(6) << x[i] << " , " << u(x[i]) << " , " << v[i] <<endl;
        }
        ofile.close();
    }
    // Adds number of points, time usded and maximum error to "run_info_file"
    ofile.open(run_info_file, ios_base::app);
    ofile << showpoint << setprecision(6) << setw(6) <<"N = " << k << " t = " << t << " eps = " << error << endl;
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
  int runs;

  if (argc<=1){
    cout << "Error: " << argv[0] << " reads number of integration points and number of runs per method (optional, default 1)" << endl;
    exit(1);
  }

  //initialising first elements in the vectors
  if (argc == 2){
    runs = 1;
  }
  else {
    runs = atof(argv[2]);
  }

  // Set boundary values
  x[0] = 0, x[n+1] = 1;
  exact[0] = exact[n+1] = 0;
  double hh = h*h; //do this outside loops to minimize FLOPS in loops

  //doing a loop for initializing
  cout << "Initializing vectors, and calculating exact vectors" << endl;
  for(int i = 1; i < n+1; i++){
    a[i] = -1, c[i] =-1, b[i] = 2; //initialising vecotors a, b, c
    x[i] = i*h; // making vector x
    b_tilde[i] = hh*f(x[i]); // making vector b_tilde
    exact[i] = u(x[i]); // finding exact solution
  }

  cout << "\nCalculating v with standard Thomas algorithm" << endl;

  // Set the first pair of file names
  string run_info_file = "run_info_standard.txt";
  // always the same, as they contain all runs
  // including previous with different amount of points

  string outfilename;
  string number = argv[1];
  outfilename = "output_"+number+".data";

  double max_error;
  double time_spent;

  a[0]=a[n+1]=c[0]=c[n+1]=-1.0;

  for(int i = 0; i < runs; i++){
    cout << "run number " << i+1 << endl;

    start = clock();
    // calculate the solution with the standard algorithm
    solver_Thomas(a,b,c,b_tilde,v,N);
    finish = clock();
    // find the error in our numerical solution
    max_error = find_max_error(exact,v,N);

    // calculate the time it took
    time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
    for(int j=1; j<n+1; j++){
      // reset arrays, as these are changed by the solver function
      b_tilde[j] = hh*f(x[j]);
      b[j] = 2;
    }
    print_file(max_error,time_spent,v,x, outfilename,n,run_info_file);
  }
  // free up space since I don't need a or c any more
  delete [] a;
  delete [] c;

  cout << "Max relative error = " << max_error << endl;
  cout <<"Time elapsed standard version = " << time_spent <<" seconds" <<endl;

  // restarts, with the faster algorithm
  cout << "\nCalculating v with specialized Thomas algorithm"<<endl;
  for(int i = 1; i <= n; i++){
    b[i] = (i+1)/(double)i; //precalculated values for b
    b_tilde[i] = hh*f(x[i]); // resets b_tilde
  }

  outfilename = "output_quick_"+number+".data";
  run_info_file = "run_info_specialized.txt";

  for(int i = 0; i<runs; i++){
    cout << "Run number " << i << endl;
    start = clock();
    // calculate the solution with the standard algorithm
    solver_Thomas_Specialized(b,b_tilde,v,N);
    finish = clock();
    // find the error in our numerical solution
    double max_error = find_max_error(exact,v,N);

    // calculate the time it took
   double time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
   for(int j=1; j<n+1; j++){
    // reset b_tilde array, as it is changed by the solver function
    b_tilde[j] = hh*f(x[j]);
    }
    print_file(max_error,time_spent,v,x, outfilename,n,run_info_file);
  }

  cout << "Max relative error = " << max_error << endl;
  cout <<"Time elapsed specialized version = " << time_spent <<" seconds" <<endl;

  // delete arrays to save space
  delete [] b;
  delete [] b_tilde;

  cout<<"\nCalculating v using the LU-decomposition"<<endl;
  // timing the LU decompositio

  // initialize matrix and vectors
  mat A = zeros<mat>(n,n);
  vec b_tilde2(n);
  vec solution(n);
  for (int i = 0; i < n - 1; i++){ //making the tridiagonal matrix
    A(i,i) = 2;
    A(i,i+1) = -1;
    A(i+1,i) =-1;
    b_tilde2[i] = hh*f(x[i+1]);
  }
  b_tilde2[n-1] = hh*f(x[n]);
  A(n-1,n-1) = 2; // setting the last corner element to -2;


  start = clock();
  solution = LUdecomposition(b_tilde2, A);
  finish = clock();

  for (int i = 0; i<n; i++){
    // transfer the numbers in the 'solution' vector to an array of floats (v)
    v[i+1] = solution[i];
  }
  run_info_file = "run_info_LU.txt";

  max_error = find_max_error(exact,v,N);
  cout << "Max relative error = " << max_error << endl;
  time_spent = ((double)(finish - start)/ CLOCKS_PER_SEC);
  cout <<"Time elapsed LU = " << time_spent << " seconds" <<endl;
  outfilename = "output_LU_"+number+".data";
  print_file(max_error,time_spent,v,x, outfilename,n,run_info_file);
  return 0;
}
