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



void writingfunc1D(int n, vec u_t, ostream& ofile)
{
ofile << setiosflags(ios::showpoint | ios::uppercase);
for (int k = 0; k <= n; k++)
  {
    ofile << setw(15) << setprecision(8) << u_t(k); // writes u for all x
  }
  ofile << endl;
}



void writingfunc2D(int n, mat u_t, ostream& ofile)
{
ofile << setiosflags(ios::showpoint | ios::uppercase);

for (int i = 0; i <= n; i++)
  {
    ofile << endl;
    for (int k = 0; k <= n; k++)
    {
    ofile << setw(15) << setprecision(8) << u_t(i,k); // writes u for all x
    }
  }
  ofile << endl;
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
}

// explicit 2D
void explicitsch2D(int n, int t_steps)
{
// setting up steplength in x = y and t
double dx = 1.0 / (n+1); //dy = dx , alpha is the same for both as well
double dt = 1.0 /(t_steps);
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

ofile << setiosflags(ios::showpoint | ios::uppercase);
ofile << n+1 << endl;

writingfunc2D(n, u_t, ofile); // to write first line at i,j = 0

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

    writingfunc2D(n, u_t, ofile);
  }

}

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
  //explicitsch1D(n,t);
  explicitsch2D(n,t);
  ofile.close();
  return 0;
}
