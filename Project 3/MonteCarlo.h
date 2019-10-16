#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>

using namespace std;

typedef double(*FunctionPointers)(double);

// Preforms MonteCarlo integration
// Takes:
// a function to integrate
// a random number function, should take a number between 0 and 1 and return a double
// The number of points to calculate
// and the dimensions of the function
double MonteCarlo(double (*)(double*), FunctionPointers*, int, int);

// Class to transform a random number between 0 and 1, to one between a and b
class uniform_transform {
  double a; // start of the uniform distribution
  double b; // stop of the uniform distribution
  double diff; // b-a, saved as variable to save time
public:
  uniform_transform(double, double);
  double transform(double);
};

#endif