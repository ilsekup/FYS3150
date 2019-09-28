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
double maxoff(mat&, mat&, int, int *, int *);
void jacobi_method(mat&, int, int, int);
void iterative(mat&, mat&, int);
mat initialize(int, double, bool);
void print_file(mat&, vec&, double, int);
double get_error(vec, vec, int);
vec analytic_eigvals_harmonic(int);
vec sort_eigenproblem(mat &, mat &, int);

#endif