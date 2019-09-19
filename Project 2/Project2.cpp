#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>

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
double maxoff(mat A, int N, int * k, int * l ){
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
mat jacobi_method(mat A, int k, int l, int N){
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
  return A;
}

int main(int argc, char *argv[]){
  int N = 3;
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
  vec eigvals = eig_sym(A);
  eigvals.print("eigvals = ");
  vec eigvals_a = analytic_eigvals(N,d,a);
  eigvals_a.print("Analytical eigenvalues = ");

  //doing the jacobimethod multiple times untill the non diagonal elements are close enough to zero
  int counter = 0;
  int k, l;
  double A_max = maxoff(A, N, &k, &l);

  while(A_max > 1e-8){
    A_max = maxoff(A, N, &k, &l);
    A = jacobi_method(A, k, l, N);
    A.print("A= ");
    counter++;
  }
  A.print("A= ");
  cout<<"number of iterations needed: "<< counter <<endl;
  return 0;
}