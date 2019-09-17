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

vec analytic_eigvals(int N,double d,double a){
  vec temp = zeros<vec>(N);
  for(int i = 0;i < N; i++){
    temp(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
  }
  return temp;
}

int main(int argc, char *argv[]){
  int N = 10;
  mat A = zeros<mat>(N,N);
  double h = 1.0/(double) N;
  double hh = h*h;
  double d = 2.0/hh;
  double a = -1/hh;

  for (int i = 0; i < N-1; i++){ //making the tridiagonal matrix
    A(i,i) = d;
    A(i,i+1) = A(i+1,i) = a;
  }
  A(N-1,N-1) = 2.0/hh;
  vec eigvals = eig_sym(A);
  eigvals.print("eigvals = ");
  vec eigvals_a = analytic_eigvals(N,d,a);
  eigvals_a.print("Analytical eigenvalues = ");

  return 0;
}

double off(mat A, int N){
  int k, l;
  double A_max = 0.;
  for(int i = 0; i<N; i++){
    for(int j=0; j<N; j++){
      if(i!=j){
        if(A_max<A(i,j)){
          A_max = A(i,j);
          k = i;
          l = j;
        }
      }
    }
  }
  return A_max,k,l;
}