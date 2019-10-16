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

// Preforms MonteCarlo integration
// Takes:
// a function to integrate
// a random number function, should take a number between 0 and 1 and return a double
// The number of points to calculate
// and the dimensions of the function
double MonteCarlo(double (*)(double* ), double* , double* , int, int);

#endif