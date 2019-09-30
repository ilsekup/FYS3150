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
  clock_t start, finish;
  double time_spent;
  int N = atof(argv[1]);
  double max = 4.0; //when we dont have a potential rho[N-1] = 1, when we do rho[N-1] =infinity
  mat A = initialize(N, max, true);
  mat R(N,N,fill::eye);
  double h = 1.0/N; //rho_N = 1 rho_0 = 0
  double hh = h*h;
  double d = 2.0/hh; //diagonal elements
  double a = -1/hh;
  //A.print("A = ");
  start = clock();
  vec eigenvalues_arma = eig_sym(A);
  finish = clock();
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
  cout << "time for amamadillo eigenvlaues: " << time_spent << " s" << endl;
  //eigenvalues_arma.print("Eigenvalues (from armadillo) = ");
  start = clock();
  vec eigenvalues_analytical = analytic_eigvals(N,d,a); //this only shows the eigenvalues with no potential
  finish = clock();
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
  cout << "Time to calculate analytically: " << time_spent << " s" << endl;
  eigenvalues_analytical.print("Analytical eigenvalues = ");
  vec eigenvalues_analytical_potential = analytic_eigvals_harmonic(N);
  start = clock();
  iterative(A,R,N);
  vec eigenvalues_Jacobi = sort_eigenproblem(A, R, N);
  eigenvalues_Jacobi.print("Eigenvalues (Jacobi) = ");
  finish = clock();
  time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
  cout << "Time spend with jacobi_method: " << time_spent << " s" << endl;
  //R.print("R = ");
  //A.print("A = ");
  double average_error = get_error(eigenvalues_analytical_potential, eigenvalues_Jacobi, N);
  cout<<"The average error for the ten first elements for N = " << N << " is: "<< average_error << endl;
  print_file(R,eigenvalues_Jacobi, average_error, N);

}
