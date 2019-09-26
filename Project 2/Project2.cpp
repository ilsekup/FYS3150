#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>

#include "Project2.h"

std::ofstream ofile;
using namespace std;
using namespace arma;

vec analytic_eigvals(int N,double d,double a){ //making a function for finding the analytical eigenvalues given in exercise
  vec temp = zeros<vec>(N);
  for(int i = 0;i < N; i++){
    temp(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
  }
  return temp;
}

//finds the largest of the non-diagonal elements in A and the position for this element
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

//prints eigenvectors (R) and eigenvalues (lambda) to files
void print_file(mat& R, vec& lambdas, int N){
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
    ofile << showpoint << setprecision(6) << setw(6) << "eigenvalues (corresponds to columns in 'eigenvectors_*_.txt')" << endl;
    ofile << "------------------------------" << endl;
    lambdas.print(ofile);
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


int main(int argc, char *argv[]){
  int N = atof(argv[1]);
  double max = 10.0; //when we dont have a potential rho[N] = 1
  mat A = initialize(N, max, true);
  mat R(N,N,fill::eye);
  double h = 1.0/N; //rho_N = 1 rho_0 = 0
  double hh = h*h;
  double d = 2.0/hh; //diagonal elements
  double a = -1/hh;
  //A.print("A = ");
  vec eigenvalues_arma;
  mat eigenvectors_arma;
  eig_sym(eigenvalues_arma, eigenvectors_arma,A);
  eigenvalues_arma.print("Eigenvalues (from armadillo) = ");
//  vec eigenvalues_analytical = analytic_eigvals(N,d,a); //this only shows the eigenvalues with no potential
//  eigenvalues_analytical.print("Analytical eigenvalues = ");
  //solve iteratively with Jacobi
  iterative(A,R,N);
  vec eigenvalues_Jacobi = get_eigvals(A,N);
  sort_eigenproblem(R,eigenvalues_Jacobi,N);
  print_file(R,eigenvalues_Jacobi,N);
}

