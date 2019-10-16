#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>
#include <stdlib.h>
#include"Gauss.h"
using namespace std;

void gauleg(double, double, double *, double *, int);
double gauss_laguerre(double*, double *, int, double);

//the function we want to integrate
double int_function(double x1, double y1, double z1, double x2, double y2, double z2){
  double alpha = 2.;
  double e1 = -2*alpha*sqrt(x1*x1+ y1*y1 + z1*z1);
  double e2 = -2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2);
  double denomenator = sqrt(pow((x1- x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2));
  //need to account for r1 - r2 = 0 because we then have division by zero.
  if(denomenator == 0){
    return 0;
  }
  else{
    return exp(e1+e2)/denomenator;
  }
  }

//calculating the sum as an approximation of the integral using weights found with legendre polynomial
double Gausslegendre(int N, double lam){
  double *x = new double [N];
  double *w = new double [N];

  gauleg(-lam, lam, x, w, N); //this function already scales for differnet integration numbers.
  double int_gauss = 0.;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      for(int k=0; k<N; k++){
        for(int l=0; l<N; l++){
          for(int m=0; m<N; m++){
            for(int n=0; n<N; n++){
              int_gauss += w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i], x[j], x[k], x[l], x[m], x[n]);

            }
          }
        }
      }
    }
  }
  cout << "Legendre method with N= "<< N << endl;
  cout <<"Integral: " <<int_gauss << endl;
  double Error = 0.19277 - int_gauss;
  cout <<"Error: " << Error << endl;
  delete [] x;
  delete [] w;
  return int_gauss;
}


//the function we want to integrate over in spherical coordinates
double int_func_sphere(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  double numerator = sin(theta1)*sin(theta2);
  double cosb = cos(theta1)*cos(theta2) + (sin(theta1)*sin(theta2)*cos(phi1-phi2));
  double denominator = sqrt(r1*r1 + r2*r2 - (2*r1*r2*cosb));
  double tol = 1e-6; //making sure the denominator is not zero or almost zero
  if(denominator < tol || r1*r1 + r2*r2 - (2*r1*r2*cosb) < 0){
    return 0.0;
  }
  return numerator/denominator;
}


//calculating the sum by finding the phi and theta weights with legendre and r weights with laguerre
double Gausslaguerre(int N){
  double *rlag = new double[N+1]; //laguerre starts from 1
  double *wlag  = new double[N+1];
  double *phi = new double[N];
  double *theta = new double[N];
  double *wphi = new double[N];
  double *wtheta = new double[N];
  double alpha = 2.0;
  gauss_laguerre(rlag, wlag, N, alpha);
  gauleg(0, 2*M_PI, phi, wphi, N);
  gauleg(0, M_PI, theta, wtheta, N);
  double int_gausslag = 0;
  for(int i=1; i<N+1; i++){
    for(int j=1; j<N+1; j++){
      for(int k=0; k<N; k++){
        for(int l=0; l<N; l++){
          for(int m=0; m<N; m++){
            for(int n=0; n<N; n++){
              int_gausslag += wlag[i]*wlag[j]*wtheta[k]*wtheta[l]*wphi[m]*wphi[n]*int_func_sphere(rlag[i], rlag[j], theta[k], theta[l], phi[m], phi[n]);
            }
          }
        }
      }
    }
  }
  //the integral of r1^2 r2^2 exp(-4(r1 + r2)) is 1/1024, add this at the end.
  int_gausslag = (1/1024.*int_gausslag);
  cout << "Laguerre method with N= "<< N << endl;
  cout <<"Integral: " <<int_gausslag << endl;
  double Errorlag = 0.19277 - int_gausslag;
  cout <<"Error: " << Errorlag << endl;
  delete [] rlag;
  delete [] wlag;
  delete [] wtheta;
  delete [] theta;
  delete [] phi;
  delete [] wphi;
  return int_gausslag;
}
