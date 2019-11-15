#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
using namespace std;

void initialize(int, double, int **, double&, double&, long, bool);
inline int periodic(int, int, int);
void metropolis(int, long&, int**, double&, double&, double*);
void writingfunc(int, int, double, double*, ostream&);

#endif