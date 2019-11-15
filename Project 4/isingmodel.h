#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>

void initialize(int, double, int **, double&, double&, long, bool);
inline int periodic(int, int, int);
void metropolis(int, long&, int**, double&, double&, double*);
void writingfunc(int, int, double, double*);
void writingfunc2(int , int , double , double*, double *);
void writingfunc3(int, int, double, double, double *, int);

#endif