#ifndef PROJECT2_H
#define PROJECT2_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

vec analytic_eigvals(int,double,double);
double maxoff(mat&, int, int *, int *);
void jacobi_method(mat&, int, int, int);
void iterative(mat&, int);
mat initialize(int);
vec get_eigvals(mat&, int);

#endif