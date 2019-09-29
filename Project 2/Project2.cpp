#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>

#include "Project2.h"

using namespace std;
using namespace arma;


//analytic eigenvalues for the buckling beam
vec analytic_eigvals(int N,double d,double a){ //making a function for finding the analytical eigenvalues given in exercise
  vec temp = zeros<vec>(N);
  for(int i = 0;i < N; i++){
    temp(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
  }
  return temp;
}

//analytic eigenvlaues for the harmonic oscillator
vec analytic_eigvals_harmonic(int N){
  vec temp = zeros<vec>(N);
  for(int i = 0; i < N; i++){
    temp(i) = 3 + 4*i;
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

//prints results to file
void print_file(mat& R, vec& lambda, double error, int N){
    std::ofstream ofile;
    string outfilename;
    string number = to_string(N);
    cout << number << endl;

    outfilename = "eigenvectors_"+number+".txt";
    ofile.open(outfilename);
    ofile << showpoint << setprecision(6) << setw(6) << "eigenvectors (coloumn format)" << endl;
    ofile << "------------------------------" << endl;
    R.print(ofile);
    ofile.close();

    outfilename = "eigenvalues_"+number+".txt";
    ofile.open(outfilename);
    ofile << "average error: " << error << endl;
    ofile << showpoint << setprecision(6) << setw(6) << "eigenvalues (corresponds to columns in 'eigenvectors_*_.txt')" << endl;
    ofile << "------------------------------" << endl;
    lambda.print(ofile);
    ofile.close();
}

//Preforms one step of the Jacobi algorithm
//this function changes the elements in the matrices A and R
void jacobi_method(mat &A, mat &R, int k, int l, int N){
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

    //Rotating the eigenvector matrix
    }
    double Rik = R(i,k);
    double Ril = R(i,l);
    R(i, k) = Rik*c - Ril*s;
    R(i, l) = Ril*c + Rik*s;
  }
}

//doing the jacobimethod multiple times untill the non diagonal elements are close enough to zero
//overwirtes matrix A
void iterative(mat &A, mat &R, int N){
  int counter = 0;
  int k, l;
  double A_max = maxoff(A, N, &k, &l);

  while(A_max > 1e-8){
    A_max = maxoff(A, N, &k, &l);
    jacobi_method(A, R, k, l, N);
    counter++;
  }
  cout<<"number of iterations needed: "<< counter <<endl;
}


//Finds the error between the analytical and computed eigenvalues for the first 20 eigenvalues
double get_error(vec analytical, vec eigenvals, int N){
  vec error = zeros<vec>(N);
    for(int i= 0; i < 20; i++){
      if(N >=20){
      error(i) = fabs((eigenvals[i] - analytical[i])/analytical[i]);
    }
    else{
      break;
  }
  }
  double average_error = accu(error)/N;
  return average_error;
}

//Returns a vector of the diagonal elements of matrix A
//Sorts a matrix R containing the eigenvectors, and a vector lambdas containing the eigenvalues
//overwirtes R and lambdas with sorted versions
vec sort_eigenproblem(mat &A, mat &R, int N){
  vec lambdas = zeros<vec>(N);
  for(int i = 0;i<N; i++){
    lambdas(i) = A(i,i);
  }
  vec lambda_new = sort(lambdas);
  mat R_new = zeros<mat>(N,N);
  ucolvec ind;
  for(int i=0; i<N; i++){
    ind = find(lambdas==lambda_new(i));
    R_new.col(i) = R.col(ind(0));
  }
  R = R_new;
  lambdas = lambda_new;
  return lambdas;
}

//Returns the tridiagonal matrix A needed to solve the boundary problem

mat initialize(int N, double max, bool potential){
  mat A = zeros<mat>(N,N);
  double *rho = new double[N];
  double *d = new double[N];
  rho[N-1] = max;
  double h = rho[N-1]/(double) N; //rho_N = max  rho_0 = 0
  double hh = h*h;
  double a = -1/hh;

  //this will only work if we have a potential, need to figure out for no potential
  for(int i = 0; i <= N-1; i++){
    if(potential == true){
      rho[i] = i*h;
    }
    else{
      rho[i] = 0;
    }
      d[i] = 2.0/hh + rho[i]*rho[i]; //diagonal elements
    }

  for (int i = 0; i < N-1; i++){ //making the tridiagonal matrix
    A(i,i) = d[i];
    A(i,i+1) = A(i+1,i) = a;
  }
  A(N-1,N-1) = d[N-1]; //setting last diagonal element, as for loop does not go this far

  //delete unneeded arrays
  delete [] rho;
  delete [] d;
  return A;
}

//the main function is a different file to avoid conflict when compiling the test functions
