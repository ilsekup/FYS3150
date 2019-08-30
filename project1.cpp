#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>

std::ofstream ofile;

int f(double x){
  return 100*exp(-10*x);
}

double u(double x){
  return 1 - (1- exp(-10))*x - exp(-10*x);
}
int main(int argc, char *argv[]){
  //initialize vectors for g, gt, d and dt and e
  int n = atof(argv[1]);
  double *v = new double[n];
  double *d = new double[n];
  double *b = new double[n];
  double *b_tilde = new double[n]; //does this have to be a n length or n+1 length ?
  double *g = new double[n];
  double *c = new double[n];
  double *a = new double[n];
  double *x = new double[n];
  double h = 1.0/(n+1);
  x[0] = 0;
  x[n+1] = 1;
  d[1] = b[1];
  g[1] = b_tilde[1];
  v[0] = v[n+1] = 0;

//doing a loop for finding the new values for vectors
  std::cout << "Making new vectors for b and b-tilde \n";
  for(int i = 1; i < n; i++){
    a[i] = -1, c[i] =-1, b[i] = 2; //initialising vecotors a, b, c
    x[i] = x[0] + i*h; // making vector x
    b_tilde[i] = h*h*f(x[i]); // making vector b_tilde
    d[i] = b[i] - (a[i-1]*c[i-1])/b[i-1]; //vector d is the new vector b
    g[i] = b_tilde[i]  - (b_tilde[i-1]*a[i-1])/b[i-1]; // vector g is the new vector b_tilde

  }
  delete [] a; // free up space since I don't need x any more
  delete [] b_tilde;
  delete [] b;
//making a decreasing for loop for finding the values for vector v
  std::cout<<"Finding the calculated values for v \n";
  for(int i = n-1; i >= 1; i--){ //loop over n-1 ? n or n+1
    v[i] = (g[i] - c[i]*v[i+1])/d[i];

  }
  delete [] c;

  //finding the error at a random point
  std::cout<<"calculating the error between calculated and expected answer \n";
  double error;
  error = u(x[50] - v[50]);
  std::cout << error << "\n";

  using namespace std;
  char *outfilename;
  outfilename = argv[2];
  ofile.open(outfilename);
  ofile << setw(15) << setprecision(8) << "relative error=" << error << endl;
  ofile.close();
}
