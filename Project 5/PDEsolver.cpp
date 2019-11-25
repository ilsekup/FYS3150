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



void writingfunc(int n,int& j, vec u_t, ostream& ofile)
{
ofile << setiosflags(ios::showpoint | ios::uppercase);
for (int k = 0; k <= n; k++)
  {
    ofile << setw(15) << setprecision(8) << u_t(k); // writes u for all x
  }
  ofile << setw(15) << setprecision(8) << j; // temp iteration, last element in row
  ofile << setw(15) << setprecision(8) << '\n' << endl;
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

for (int j = 1; j < t_steps; j++) // iterating over temperatures
  {
    for (int i = 1; i < n; i++) // iterating over x-position
    {
    u_t(i) = alpha*( u_xx(i-1) + u_xx(i+1) ) + (1 - 2*alpha) * u_xx(i);
    u_xx(i) = u_t(i);
    }
    writingfunc(n,j, u_t, ofile);
  }

}

int main(int argc, char* argv[])
{
  char *outfilename;
  outfilename = argv[1];
  ofstream ofile;
  ofile.open("explicit.txt");
  int n, t; // number of steps in x and t respectively.

  cout << "n = " << endl;
  cin >> n;
  cout << "t = " << endl;
  cin >> t;
  explicitsch(n,t);

  ofile.close();
  return 0;
}
