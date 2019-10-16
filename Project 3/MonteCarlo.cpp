#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include "MonteCarlo.h"
#include <random>
#include <functional>

using namespace std;
using namespace std::placeholders;

double MonteCarlo(std::function<double(double*)> func, double* a, double* b, int N, int d){
    uniform_real_distribution<double> distribution(0.0,1.0);
    ranlux48 generator;
    generator.seed(time(NULL));
    auto random = bind(distribution,generator);

    // Define arrays for a-b and the random values for the position vector r
    double* diff = new double[d];
    double* r = new double[d];
    double jacobian = 1.0;

    for(int i=0;i<d;i++){
      diff[i] = b[i] - a[i];
      jacobian *= diff[i];
    }

    double sum = 0;

    for(int i=0; i<N;i++){
        for(int j=0;j<d;j++){
            r[j] = a[j] + random()*diff[j];
        }
        sum += func(r);
    }
    delete [] r;
    delete [] diff;

    sum *= jacobian;
    return sum / (double) N;
}