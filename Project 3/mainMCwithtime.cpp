#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include "MonteCarloMP.h"
#include <functional>

using namespace std;
using namespace std::placeholders;

// Function without importance sampling
double f(double* );

// y(x)
void transform_r(double*, double);

// Short class for function with imortance sampling (f(x(y))/p(x(y))), for varying alpha
class F{
  double alpha;
  double sqr_alpha;
public:
  F(double _alpha){
    alpha = _alpha;
    sqr_alpha = alpha*alpha;
  }
  double f_on_p(double* );
};
double F::f_on_p(double* x){
  transform_r(x,alpha);
  double r1sqr, r2sqr;
  r1sqr = x[0]*x[0];
  r2sqr = x[3]*x[3];
  double exp_on_p = exp((alpha-4)*(x[0]+x[3]))/sqr_alpha;
  double cosbeta = cos(x[1])*cos(x[4]) + sin(x[1])*sin(x[4])*cos(x[2]-x[5]);
  return exp_on_p*r1sqr*r2sqr*sin(x[1])*sin(x[4])/(sqrt(r1sqr + r2sqr - 2*x[0]*x[3]*cosbeta));
}

// Function to integrate when alpha = 4, less computationally heavy
// (for pard d and e)
double f2(double*);

int main(int argc, char *argv[]){

    clock_t start, finish;
    double time_spent;
    int N = atof(argv[1]);
    double exact = 5*M_PI*M_PI/256;

    double inf = 2.0;

    double* begin = new double[6];
    double* stop = new double[6];

    for(int i=0;i<6;i++){
      begin[i] = -inf;
      stop[i] = inf;
    }
    array<double, 2> MC_estimate;
    array<double, 2> MC_importance_sampling;
    start = clock();
    MC_estimate = MonteCarloMP(f,begin,stop,N,6);
    finish = clock();
    time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
    cout << "Time spent MC = "<< time_spent << " seconds"<<endl;
    for(int i=0;i<6;i++){
      begin[i] = 0;
    }
    stop[0] = stop[3] = 1;
    stop[1] = stop[4] = M_PI;
    stop[2] = stop[5] = 2*M_PI;


    // How to use the function, if one want alpha != 4:
    //double alpha = 2.0;
    //F instance(alpha);
    //auto int_func = std::bind(&F::f_on_p, instance, _1);
    start = clock();
    MC_importance_sampling = MonteCarloMP(f2,begin,stop,N,6);
    finish = clock();
    time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
    cout << "Time spent MC (IS)  = "<< time_spent << " seconds"<< endl;
    cout << "Monte Carlo estimate                        =  " << MC_estimate[0] << " sigma = " << MC_estimate[1] << endl;
    cout << "Monte Carlo estimate (importance sampling)  =  " << MC_importance_sampling[0] << " sigma = " << MC_importance_sampling[1] << endl;
    cout << "Exact solution                              =  " << exact << endl;
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

double f2(double* x){
  transform_r(x,4.0);
  double r1sqr, r2sqr;
  r1sqr = x[0]*x[0];
  r2sqr = x[3]*x[3];
  double cosbeta = cos(x[1])*cos(x[4]) + sin(x[1])*sin(x[4])*cos(x[2]-x[5]);
  return 0.0625*r1sqr*r2sqr*sin(x[1])*sin(x[4])/(sqrt(r1sqr + r2sqr - 2*x[0]*x[3]*cosbeta));
}

void transform_r(double* r, double alpha){
  r[0] = -log(1-r[0])/alpha;
  r[3] = -log(1-r[3])/alpha;
}
