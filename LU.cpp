#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

int main (int argc, char* argv[])
{

int n = atof(argv[1]); // matrix/vector size
mat A = zeros<mat>(n,n);

for (int i = 0; i < n - 1; i++){ //making the tridiagonal matrix
  A(i,i) = -2;
  A(i,i+1) = 1;
  A(i+1,i) =1;
}

A(n-1,n-1) = -2; // setting the last corner element to -2;

A.print("A");

mat L = zeros<mat>(n,n);
mat U = zeros<mat>(n,n);

lu(L,U,A); // LU-decomposition function from armadillo

L.print(" L= ");
U.print(" U= ");
(A - L*U).print("Test"); // should print a matrix filled with zeros

return 0;
}
