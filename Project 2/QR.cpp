//Lanczos'-methods
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

//Returns a vector of the diagonal elements of matrix A
vec get_eigvals(mat &A, int N){
   vec lambdas = zeros<vec>(N);
   for(int i = 0;i<N; i++){
     lambdas(i) = A(i,i);
   }
   //eigvals = sort(eigvals);
   return lambdas;
}

//Sorts a matrix R containing the eigenvectors, and a vector lambdas containing the eigenvalues
//overwirtes R and lambdas with sorted versions
void sort_eigenproblem(mat &R, vec &lambda, int N){
  vec lambda_new = sort(lambda);
  mat R_new = zeros<mat>(N,N);
  ucolvec ind;
  for(int i=0; i<N; i++){
    ind = find(lambda==lambda_new(i));
    R_new.col(i) = R.col(ind(0));
  }
  R = R_new;
  lambda = lambda_new;
}

// K is the number of times we run the QR-factorization on mat A
// QR-factorization of tridiagonal matrix A, k-times
void QR(mat &A, int k){
  mat Q,R;
  for(int i = 0; i < k; i++){
    qr_econ(Q,R,A); // returns Q and R matrix that make A
    A = R*Q;  // Q^-1 * A * Q = A_1 = R*Q, rinse and repeat for A_k
  }
}

//Proformes m step of the Lanczos method to find eigenvalues
//Returns a tridiagonal matrix T
mat Lanczos(mat &A,int N, int m){
  mat Q = randu<mat>(N,m);
  mat T = zeros<mat>(m,m);
  vec alpha = zeros<vec>(m);
  vec r = zeros<vec>(N);
  int k = 0;
  vec beta = ones<vec>(m);

  //Normalizing q_0
  Q.col(0) /= norm(Q.col(0),2);
  r = Q.col(0);
  alpha(0) = dot(Q.col(0).t(),A*Q.col(0));
  //Do m steps
  while(beta(k)!=0 && k<m-1){
    Q.col(k+1) = r/beta(k);
    k+=1;
    alpha(k) = dot(Q.col(k).t(),A*Q.col(k));
    r = (A - alpha(k)*eye<mat>(N,N))*Q.col(k) - beta(k-1)*Q.col(k-1);
    beta(k) = norm(r,2);
  }
  //Construct the matrix T
  for(int i = 1;i<m;i++){
    T(i-1,i-1) = alpha(i-1);
    T(i,i-1) = T(i-1,i) = beta(i-1);
  }
  T(m-1,m-1) = alpha(m-1);
  if(k == m-1){
    cout << "Max iterations reached" << endl;
  }
  else{
    cout << "Norm of r was 0, aborted algorithm" << k <<endl;
  }
  return T;
}

// borrowed some main code to make tri-diag. matrix and test QR-algo.
int main(int argc, char *argv[]){
  int N = atof(argv[1]);
  int m = 10;
  if(argc<3){
    m = N;
  }
  else{
    m = atof(argv[2]);
  }
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

  //Preforme Lanczos on A
  mat T = Lanczos(A,N,m);
  //Find the eigenvalues of T with the QR algorithm
  QR(T,20);
  vec temp = analytic_eigvals(N, d, a);
  vec eigvals_Lanczos = get_eigvals(T,m);
  mat standin = randu<mat>(m,m);//dummy matrix passed to sort to get eigenvalues
  sort_eigenproblem(standin,eigvals_Lanczos,m);
  cout << "Analytical Eigenvalues (" << m <<" first) = " << endl;
  for(int i=0;i<m;i++){
    cout << "  " << temp(i) << endl;
  }
  eigvals_Lanczos.print("Eigenvalues Lanczos = ");
  return 0;
}
