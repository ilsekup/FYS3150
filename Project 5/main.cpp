#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "PDEsolver.h"
#include "lib.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  //setup for writing in to file
  char *outfilename;
  outfilename = argv[1];

  //choosing n ad t steps
  int n, t; // number of steps in x and t respectively.
  cout << "n = " << endl;
  cin >> n;
  cout << "t = " << endl;
  cin >> t;

  //calling function which also does the writing into file
  explicitsch1D(n,t);
  implicit(n, t);
  explicitsch2D(n,t);
  return 0;
}
