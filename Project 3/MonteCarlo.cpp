#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>


using namespace std;

// Preforms MonteCarlo integration
// Takes:
// a function to integrate
// a random number function, should take a number between 0 and 1 and return a double
// The number of points to calculate
// and the dimensions of the function

double MonteCarlo(double (*)(double*), double (*)(), int, int);

// The function to integrate
double func(double*);

// Placeholder function that generates random uniform numbers
// Should be rewritten, maybe as a class
double uniform();

//set engine
ranlux48 engine;
double inf = 5.0;

int main(int argc, char *argv[]){
    int N = atof(argv[1]);
    double exact = 5*M_PI*M_PI/(256);
    
    engine.seed(time(NULL));

    double MC_estimate = MonteCarlo(func,uniform,N,6);
    MC_estimate *= pow((2*inf),6);

    cout << "Monte Carlo estimate    =  " << MC_estimate << endl;
    cout << "Excat solution estimate =  " << exact << endl;
    return 0;
}

double MonteCarlo(double (*func)(double*), double (*random)(), int N, int d){
    double* r = new double[d];
    double sum = 0;

    for(int i=0; i<N;i++){
        for(int j=0;j<d;j++){
            r[j] = random();
        }
        sum += func(r);
    }
    return sum/ (double) N;
}

double func(double* r){
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

double uniform(){
    uniform_real_distribution<double> distribution(-inf,inf);
    return distribution(engine);
}
