#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include "MonteCarlo.h"
#include <functional>

using namespace std;
using namespace std::placeholders;

// Function without importance sampling
double f(double* );

// y(x)
void transform_r(double*, double);

// Prints the number of points, estimate of the integral, time spent and the standard deviation of the estimate
// Takes the method used as a string, and saves result as estimates_(name).txt
void print_file(int, double, double, double, string);

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
    using namespace std::chrono;
    // Exact value
    double exact = 5*M_PI*M_PI/256;
    // set estimate for infinity
    double inf = 2.0;
    // Set number of times to test each Monte Carlo method
    int num_tries = 100;
    int N = atof(argv[1]);

    double* start = new double[6];
    double* stop = new double[6];

    for(int i=0;i<6;i++){
      start[i] = -inf;
      stop[i] = inf;
    }
    array<double, 2> MC_estimate;
    high_resolution_clock::time_point t1, t2;
    duration<double, ratio<1,1>> t;

    string method = "MonteCarlo";
    unsigned int seed = 1;

    for(int i=0;i<num_tries;i++){
      cout << "Starting brute force run " << i << endl;
      t1 = high_resolution_clock::now();
      MC_estimate = MonteCarlo(f,start,stop,N,6,seed);
      t2 = high_resolution_clock::now();
      t = t2 - t1;
      seed += 10;
      print_file(N,MC_estimate[0],MC_estimate[1],t.count(),method);
    }


    for(int i=0;i<6;i++){
      start[i] = 0;
    }
    stop[0] = stop[3] = 1;
    stop[1] = stop[4] = M_PI;
    stop[2] = stop[5] = 2*M_PI;


    // How to use the function, if one want alpha != 4:
    //double alpha = 2.0;
    //F instance(alpha);
    //auto int_func = std::bind(&F::f_on_p, instance, _1)
    method = "MonteCarlo_imps"; // short for MonteCarlo importance sampling
    for(int i=0;i<num_tries;i++){
      cout << "Starting importance sampling run " << i << endl;
      t1 = high_resolution_clock::now();
      MC_estimate = MonteCarlo(f2,start,stop,N,6,seed);
      t2 = high_resolution_clock::now();
      t = t2 - t1;
      seed += 10;
      print_file(N,MC_estimate[0],MC_estimate[1],t.count(),method);
    }
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
  double r1sqr, r2sqr, sintheta1_sintheta2;
  r1sqr = x[0]*x[0];
  r2sqr = x[3]*x[3];
  sintheta1_sintheta2 = sin(x[1])*sin(x[4]);
  double cosbeta = cos(x[1])*cos(x[4]) + sintheta1_sintheta2*cos(x[2]-x[5]);
  return 0.0625*r1sqr*r2sqr*sintheta1_sintheta2/(sqrt(r1sqr + r2sqr - 2*x[0]*x[3]*cosbeta));
}

void transform_r(double* r, double alpha){
  r[0] = -log(1-r[0])/alpha;
  r[3] = -log(1-r[3])/alpha;
}

void print_file(int N, double I, double sigma, double t, string method){
    std::ofstream ofile;
    string outfilename;

    outfilename = "estimates_"+method+".txt";
    ofile.open(outfilename, ios_base::app);
    ofile << showpoint << setprecision(6) << setw(6) <<"N = " << N << " I = " << I << " sigma = " << sigma << " t = " << t << endl;
    ofile.close();
}