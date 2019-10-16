#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include "MonteCarlo.h"

using namespace std;

// Function without importance sampling
double f(double* );

// y(x)
void transform_r(double*, double);

// Short class for function with imortance sampling (f(x(y))/p(x(y))), as lambda may vary
class F{
  double lambda;
  double sqr_lambda;
public:
  F(double _lambda){
    lambda = _lambda;
    sqr_lambda = lambda*lambda;
  }
  double f_on_p(double* );
};
double F::f_on_p(double* x){
  transform_r(x,lambda);
  double _f = f(x);
  return sqr_lambda*exp(-lambda*(x[0] + x[3]));
}


int main(int argc, char *argv[]){
    int N = atof(argv[1]);
    double exact = 5*M_PI*M_PI/256;

    double inf = 5;
    //uniform_transform infinity(-inf,inf);

    double* start = new double[6];
    double* stop = new double[6];

    for(int i=0;i<6;i++){
      start[i] = -inf;
      stop[i] = inf;
    }
    double MC_estimate = 0;
    double MC_importance_sampling = 0;

    MC_estimate = MonteCarlo(f,start,stop,N,6);

    for(int i=0;i<6;i++){
      start[i] = 0;
    }
    stop[0] = stop[3] = 1;
    stop[1] = stop[4] = M_PI;
    stop[2] = stop[5] = 2*M_PI;

    double lambda = 2.0;

    F Func(lambda);

    MC_importance_sampling = MonteCarlo(Func.f_on_p,start,stop,N,6);

    cout << "Monte Carlo estimate                        =  " << MC_estimate << endl;
    cout << "Monte Carlo estimate (impottance sampling)  =  " << MC_importance_sampling << endl;
    cout << "Excat solution estimate                     =  " << exact << endl;
    return 0;
}

double f(double* r){
    double normr1, normr2, diff;

    normr1 = normr2 = diff = 0;

    for(int i=0;i<3;i++){
        normr1 += r[2*i]*r[2*i];
        normr2 += r[2*i+1]*r[2*i+1];
        diff += (r[2*i+1] - r[2*i])*(r[2*i+1] - r[2*i]);
    }
    normr1 = sqrt(normr1);
    normr2 = sqrt(normr2);
    diff = sqrt(diff);
    return exp(-4*(normr1 + normr2)) / diff;
}

void transform_r(double* r, double lamb){
  r[0] = -log(1-r[0])/lamb;
  r[3] = -log(1-r[3])/lamb;
}