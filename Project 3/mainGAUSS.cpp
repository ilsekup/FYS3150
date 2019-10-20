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
   double integral, weight;
  //  tie(integral, weight) = Gausslegendre(N, lam)
  //  tie(integrallag, weightlag)= Gausslaguerre(N);
   clock_t start, finish;
   double time_spent;
   for(int i = 1; i <= 3; i++){
     start = clock();
     double int_gauss, weight;
     tie(int_gauss, weight) = Gausslegendre(N*i, lam);
     finish = clock();
     time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
     cout << "Legendre method with N= "<< N*i << endl;
     cout <<"Integral: " <<int_gauss << endl;
     double Error = 0.19277 - int_gauss;
     cout <<"Error: " << Error << endl;
     cout << "Time spent for legendre with " << N*i << " iterations is " << time_spent << endl;
     start = clock();
     double int_gausslag, weightlag;
     tie(int_gausslag, weightlag) = Gausslaguerre(N*i);
     finish = clock();
     time_spent = ( (double)(finish - start)/ CLOCKS_PER_SEC );
     cout << "Laguerre method with N= "<< N*i << endl;
     cout <<"Integral: " <<int_gausslag << endl;
     double Errorlag = abs(0.19277 - int_gausslag);
     cout <<"Error: " << Errorlag << endl;
     cout << "Time spent for laguerre with " << N*i << " iterations is " << time_spent << endl;
    }

}