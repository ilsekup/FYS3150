#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <functional>

using namespace std;
using namespace std::placeholders;

// Preforms MonteCarlo integration
// Takes:
// a function to integrate
// a random number function, should take a number between 0 and 1 and return a double
// The number of points to calculate
// and the dimensions of the function
array<double, 2> MonteCarlo(std::function<double(double*)>, double* , double* , int, int, unsigned int);

#endif