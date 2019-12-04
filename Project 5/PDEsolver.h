#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include<armadillo>

using namespace std;
using namespace arma;

void writingfunc1D(int, vec, ostream&);
void writingfunc2D(int,mat, ostream&);
void explicitsch1D(int,double, int);
void explicitsch2D(int,double, int);
double *implicit(int,double, double);
void solver_Thomas(double*, double *, double *, double *, double *, int);
void CN(int , double , int , ostream&);
int implicit2D(int, double, int);
int jacobi(int, double, mat &, double, int);

#endif
