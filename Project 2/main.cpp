#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>
#include "Project2.h"
#include "Project2.cpp"


int main(int argc, char *argv[]){
  using namespace std;
  using namespace arma;
  int N = atof(argv[1]);
  double max = 1.0; //when we dont have a potential rho[N-1] = 1, when we do rho[N-1] =infinity
  mat A = initialize(N, max, false);
  mat R(N,N,fill::eye);
  double h = 1.0/N; //rho_N = 1 rho_0 = 0
  double hh = h*h;
  double d = 2.0/hh; //diagonal elements
  double a = -1/hh;
  //A.print("A = ");
  vec eigenvalues_arma = eig_sym(A);
  //eigenvalues_arma.print("Eigenvalues (from armadillo) = ");
  vec eigenvalues_analytical = analytic_eigvals(N,d,a); //this only shows the eigenvalues with no potential
  //eigenvalues_analytical.print("Analytical eigenvalues = ");
  vec eigenvalues_analytical_potential = analytic_eigvals_harmonic(N);
  iterative(A,R,N);
  vec eigenvalues_Jacobi = get_eigvals(A,N);
  eigenvalues_Jacobi.print("Eigenvalues (Jacobi) = ");
  //R.print("R = ");
  //A.print("A = ");
  double average_error = get_error(eigenvalues_analytical, eigenvalues_Jacobi, N);
  cout<<"The average error for the ten first elements for N = " << N << " is: "<< average_error << endl;
  print_file(R,eigenvalues_Jacobi, average_error, N);

}
