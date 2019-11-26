#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "lib.h"
#include <armadillo>


using namespace std;
using namespace arma;

ofstream ofile("explicit.txt", ios::out);
ofstream ofile2("implicit.txt", ios::out);



void writingfunc(int n, vec u_t, ostream& ofile)
{
ofile << setiosflags(ios::showpoint | ios::uppercase);
for (int k = 0; k <= n; k++)
  {
    ofile << setw(15) << setprecision(8) << u_t(k ); // writes u for all x
  }
  ofile << endl;
}


//explicit
void explicitsch(int n, int t_steps)
{
// setting up steplength in x and t
double dx = 1.0 / (n+1);
double dt = 1.0 /(t_steps);
double alpha = dt / ( dx * dx);
if( alpha < 0.5) // stability constraint
  {
    cout << "Warning alpha is bigger than 0.5, reconsider dt and dx" << alpha << endl;
  }

// setting up vectors and initial/bouandry conditions
vec u_xx = zeros<vec>(n+1);
vec u_t = zeros<vec>(n+1);
u_xx(0) = u_t(0) = 0.0;
u_xx(n) = u_t(n) = 1.0;

writingfunc(n, u_t, ofile); // to write first line at i,j = 0
for (int j = 1; j < t_steps + 1; j++) // iterating over temperatures
  {
    for (int i = 1; i < n; i++) // iterating over x-position
    {
    u_t(i) = alpha*( u_xx(i-1) + u_xx(i+1) ) + (1 - 2*alpha) * u_xx(i);
    u_xx(i) = u_t(i);
    }
    writingfunc(n, u_t, ofile);
  }

}
void writingfunc2(int n, double *u, ostream& ofile){
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for (int k = 0; k <= n; k++)
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

void implicit(int n, double t_steps){
  double dx = 1/double(n);
  double dt = 0.5*dx*dx;
  double alpha = 0.5; //dt/(dx*dx);
  double *vold = new double[n+2];
  double *vnew = new double[n+2];
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
}

// void CN(int n){
//
// }

int main(int argc, char* argv[])
{
  //setup for writing in to file
  char *outfilename;
  outfilename = argv[1];
  ofstream ofile;
  ofile.open("explicit.txt");

  //choosing n ad t steps
  int n, t; // number of steps in x and t respectively.
  cout << "n = " << endl;
  cin >> n;
  cout << "t = " << endl;
  cin >> t;

  //calling function which also does the writing into file
  explicitsch(n,t);
  implicit(n, t);
  ofile.close();
  return 0;
}
