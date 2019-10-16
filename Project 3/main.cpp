#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <string>
#include <armadillo>
#include "Gauss.h"
#include "Gauss.cpp"

int main(){
   int N;
   double lam;
   cout << "Enter interger for N: " << endl;
   cin >> N;
   cout << "Enter double for lambda: " << endl;
   cin >> lam;
   for(int i = 1; i < 3; i++){
      double integral = Gausslegendre(N*i, lam);
      double integrallag = Gausslaguerre(N*i);
}

}