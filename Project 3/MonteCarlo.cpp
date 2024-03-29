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

array<double, 2> MonteCarlo(std::function<double(double*)> func, double* a, double* b, int N, int d, unsigned int s){
    uniform_real_distribution<double> distribution(0.0,1.0);
    mt19937_64 generator;
    generator.seed(s);
    auto random = bind(distribution,generator);

    // Define arrays for a-b and the random values for the position vector r
    double* diff = new double[d];
    double* r = new double[d];
    double jacobian = 1.0;

    // Calculate the jacobian
    for(int i=0;i<d;i++){
      diff[i] = b[i] - a[i];
      jacobian *= diff[i];
    }

    double sum = 0, sum_sqr = 0;
    double fr;
    // Find the sum f(x_i), and the sum of f(x_i)^2
    for(int i=0; i<N;i++){
        for(int j=0;j<d;j++){
            r[j] = a[j] + random()*diff[j];
        }
        fr = func(r);
        sum += fr;
        sum_sqr += fr*fr;
    }
    // Delete arrays
    delete [] r;
    delete [] diff;
    sum_sqr /= (double) N;
    sum /= (double) N;

    double sigma = jacobian*sqrt((sum_sqr - sum*sum) / (double) N);

    sum *= jacobian;
    //sum = sum / (double) N;

    return {sum,sigma};
}