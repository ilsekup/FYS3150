#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>

//#include "Project2.h"

using namespace std;
using namespace arma;

vec analytic_eigvals(int N,double d,double a){ //making a function for finding the analytical eigenvalues given in exercise
  vec temp = zeros<vec>(N);
  for(int i = 0;i < N; i++){
    temp(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
  }
  return temp;
}

//finds the largest of the non-diagonal elements and the position for this
double maxoff(mat &A, int N, int * k, int * l ){
  double A_max = 0.;
  for(int i = 0; i<N; i++){
    for(int j=0; j<N; j++){
      if(i!=j){
        if(A_max<fabs(A(i,j))){
          A_max = fabs(A(i,j));
          *k = i;
          *l = j;
        }
      }
    }
  }
  return A_max;
}

//this function changes the elements in the matrix
void jacobi_method(mat &A, int k, int l, int N){
  double tau = (A(l,l)- A(k, k))/(2*A(k,l));
  double t;
  if(tau >= 0){
    t = 1.0/(tau + sqrt(1.0 + tau*tau));
  }
  else{
    t = - 1.0/(-tau + sqrt(1.0 + tau*tau));
  }
  double c = 1/sqrt(1 + t*t);
  double s = c*t;
  double Akk = A(k, k); //so that we can use it further
  double All = A(l,l);
  A(k, k) = Akk*c*c - 2*A(k, l)*c*s + All*s*s;
  A(l, l) = All*c*c + 2*A(k, l)*c*s + Akk*s*s;
  A(k, l) = 0;
  A(l, k) = 0;
  for(int i = 0; i < N; i++){
    if(i != k && i != l){
      double Aik = A(i, k);
      double Ail = A(i, l);
      A(i, k) = Aik*c - Ail*s;
      A(i, l) = Ail*c + Aik*s;
      A(l, i) = A(i, l);
      A(k, i) = A(i, k);
    }
  }
}

//doing the jacobimethod multiple times untill the non diagonal elements are close enough to zero
//overwirtes matrix A
void iterative(mat &A, int N){
  int counter = 0;
  int k, l;
  double A_max = maxoff(A, N, &k, &l);

  while(A_max > 1e-8){
    A_max = maxoff(A, N, &k, &l);
    jacobi_method(A, k, l, N);
    counter++;
  }
  A.print("A = ");
  cout<<"number of iterations needed: "<< counter <<endl;
}

//Returns a sorted vector of the diagonal elements
vec get_eigvals(mat &A, int N){
   vec eigvals = zeros<vec>(N);
   for(int i = 0;i<N; i++){
     eigvals(i) = A(i,i);
   }
   eigvals = sort(eigvals);
   eigvals.print("Eigenvalues (Jacobi) = ");
   return eigvals;
}

mat initialize(int N){
  mat A = zeros<mat>(N,N);
  double h = 1.0/(double) N; //rho_N = 1 rho_0 = 0
  double hh = h*h;
  double d = 2.0/hh; //diagonal elements
  double a = -1/hh;

  for (int i = 0; i < N-1; i++){ //making the tridiagonal matrix
    A(i,i) = d;
    A(i,i+1) = A(i+1,i) = a;
  }
  A(N-1,N-1) = 2.0/hh; //setting last diagonal element, as for loop does not go this far
  return A;
}

int main(int argc, char *argv[]){
  int N = atof(argv[1]);
  mat A = initialize(N);
  double h = 1.0/N; //rho_N = 1 rho_0 = 0
  double hh = h*h;
  double d = 2.0/hh; //diagonal elements
  double a = -1/hh;

  vec eigenvalues_arma = eig_sym(A);
  eigenvalues_arma.print("Eigenvalues (from armadillo) = ");

  vec eigenvalues_analytical = analytic_eigvals(N,d,a);
  eigenvalues_analytical.print("Analytical eigenvalues = ");

  iterative(A, N);
  vec eigenvalues = get_eigvals(A,N);
}