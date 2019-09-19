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
//this function changes the elements in the matrix
void jacobi_method(mat A, int k, int l){
  double tau = (A(l,l)- A(k, k))/2*A(k,l);
  double t;
  if(tau > 0){
    t = 1.0/(tau + sqrt(1.0 + tau*tau));
  }
  else{
    t = - 1.0/(-tau + sqrt(1.0 + tau*tau));
  }
  double c = 1/sqrt(1 + t*t);
  double s = t*c;
  double Akk = A(k, k); //so that we can use it further
  double All = A(l,l),
  A(k, k) = Akk*c*c - 2*A(k, l)*c*s + All*s*s;
  A(l, l) = All*c*c - 2*A(k, l)*c*s + Akk*s*s;
  A(k, l) = 0;
  A(l, k) = 0;
  if(i != k && i != l){
    double Aik, Ail
    A(i, k) = Aik
    A(i, l) = Ail
    A(i, k) = Aik*c - Ail*s
    A(i, l) = Ail*c + Aik*s
  }
  return;
}

int main(int argc, char *argv[]){
  int N = 100;
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
  double A_max, k, l = off(A, N);
  cout << "max= " << k << endl;

  while(A_max > 1e-8){
    double A_max, k, l = off(A, N);
    cout<<"k= " << k << endl;
    jacobi_method(A, k, l);
    counter++;
  }
  cout<<"number of iterations needed: "<< counter <<endl;
  return 0;
}