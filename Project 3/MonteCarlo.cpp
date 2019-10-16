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

// The function to integrate
//double func(double*);

uniform_transform::uniform_transform(double a1, double b1){
  a = a1;
  b = b1;
  diff = b-a;
}

double uniform_transform::transform(double x){
  return a + x*diff;
}

double MonteCarlo(double (*func)(double*), FunctionPointers* transforms, int N, int d){
    //set engine
    uniform_real_distribution<double> distribution(0,1);
    ranlux48 generator;
    generator.seed(time(NULL));
    auto random = bind(distribution,generator);

    double* r = new double[d];
    double sum = 0;

    for(int i=0; i<N;i++){
        for(int j=0;j<d;j++){
            r[j] = (transforms[j])(random());
        }
        sum += func(r);
    }
    delete [] r;
    return sum / (double) N;
}

/*
double uniform(){
    uniform_real_distribution<double> distribution(-inf,inf);
    return distribution(engine);
}
*/