#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"

std::ofstream ofile;

double f(double x){
  return 100*exp(-10*x);
}

double u(double x){
  return 1 - (1- exp(-10))*x - exp(-10*x);
}
int main(int argc, char *argv[]){
  clock_t start, finish;
  //initialize vectors for a, b, c, b_tilde, x and v
  int n = atof(argv[1]);//number of points we calculate on
  int N = n+2;//actual number of elements (includes x_0 and x_n+1)
  double *v = new double[N];
  double *b = new double[N];
  double *b_tilde = new double[N];
  double *c = new double[N];
  double *a = new double[N];
  double *x = new double[N];
  double h = 1.0/(n+1);
  x[0] = 0, x[1] = h, x[n+1] = 1; //initialising first elements in the vectors
  a[1] = -1., b[1] = 2., c[1] = -1;
  double hh = h*h; //do this outside loop to minimize flops in loop

//doing a loop for finding the new values for vectors
  std::cout << "Starting step 1, overwrites b and b_tilde" << std::endl;
  start = clock();
  for(int i = 2; i < n+1; i++){
    a[i] = -1, c[i] =-1, b[i] = 2; //initialising vecotors a, b, c
    x[i] = i*h; // making vector x
    b_tilde[i] = hh*f(x[i]); // making vector b_tilde
    b[i] -= c[i-1]*a[i-1]/b[i-1];
    b_tilde[i] -= b_tilde[i-1]*a[i-1]/b[i-1];

  }
  delete [] a; // free up space since I don't need a any more
//making a decreasing for loop for finding the values for vector v
  std::cout<<"Finding the calculated values for v \n";
  for(int i = n-1; i>0; i--){
    b_tilde[i] -= b_tilde[i+1]*c[i]/b[i+1];
    v[i] = b_tilde[i]/b[i];
  }
  v[n] = b_tilde[n]/b[n];
  finish = clock();
  delete [] c;

  //summing the error at all points
  using namespace std;
  std::cout<<"Calculating the error between calculated and expected answer \n";

  double rel_error;
  double max_error = log10(abs((v[1] - u(x[1]))/u(x[1])));
  for (int i=1; i<n; i++){
    rel_error = log10(abs((v[i] - u(x[i]))/u(x[i]))); //logarithm of relative error
    if (rel_error>max_error){
        max_error = rel_error;
    }
    }
  cout << "Max relative error = " << max_error << endl;
  cout <<"Time elapsed =  " << ( (double)(finish - start)/ CLOCKS_PER_SEC ) <<" seconds" <<endl;
  char *outfilename;
  outfilename = argv[2];
  ofile.open(outfilename);
  ofile << "Max error = " << max_error << endl;
  ofile << "Time used = " << ( (double)(finish - start)/ CLOCKS_PER_SEC ) << endl;
  ofile << "   x  ,   u   ,   v" << endl;
  ofile << "------------------------------" << endl;
  for(int i = 0; i < n+1; i++){
    ofile << showpoint << setprecision(6) << setw(6) << x[i] << " , " << u(x[i]) << " , " << v[i] << endl;
  }
  ofile.close();
  return 0;
}
