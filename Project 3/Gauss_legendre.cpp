#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>
#include <stdlib.h>
// #define EPS 3.0e-14
// #define MAXIT 10
// #define   ZERO       1.0E-10
//#include <mpi.h>
using namespace std;
// Integral is from -inf to inf of  e^-2a(r1 + r2)
// Where r1 = sqrt(x_1^2 + y1^2 + z1^2 )
void gauleg(double, double, double *, double *, int);

double int_function(double x1, double y1, double z1, double x2, double y2, double z2){
  double alpha = 2.;
  double e1 = -2*alpha*sqrt(x1*x1+ y1*y1 + z1*z1);
  double e2 = -2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2);
  double denomenator = sqrt(pow((x1- x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2));
  return exp(e1+e2)/denomenator;
  }

int main(){
  int N = 20;
  double lam = 2;
  double *x = new double [N];
  double *w = new double [N];

  gauleg(-lam, lam, x, w, N);
  double int_gauss = 0.;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      for(int k=0; k<N; k++){
        for(int l=0; l<N; l++){
          for(int m=0; m<N; m++){
            for(int n=0; n<N; n++){
              int_gauss += w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(lam*x[i], lam*x[j], lam*x[k], lam*x[l], lam*x[m], lam*x[n]);

            }
          }
        }
      }
    }
  }
  int_gauss = lam*int_gauss;
  cout << int_gauss << endl;
}
