#ifndef GAUSS_H
#define GAUSS_H

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
void gauleg(double, double, double *, double *, int);
double int_function(double, double, double, double, double , double );
double Gausslegendre(int, double);
double gauss_laguerre(double *, double *, double, int);
double Gausslaguerre(int N);
double int_func_sphere(double, double, double, double, double, double);



#endif