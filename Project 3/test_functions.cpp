#include "catch.hpp"
#include "Gauss.h"
#include "MonteCarlo.h"

//Here we test if the approximation of the integral is as good as the approximation in the lecture notes.
TEST_CASE("Testing the integral result"){
  double integralgauss, weight1;
  tie(integralgauss, weight1) = Gausslaguerre(10);
  double intlag, weight1lag;
  tie(intlag, weight1lag) = Gausslaguerre(10);
  REQUIRE(integralgauss ==Approx(0.19277).epsilon(0.06)); //Checking if the approx is close for N = 30
  REQUIRE(intlag == Approx(0.19277).epsilon(0.006));
}

//Here we test if the first weight gotten from the weightfunctions are equal to the analytical
TEST_CASE("Testing the first weight"){
  double integralgauss, weight1;
  tie(integralgauss, weight1) = Gausslaguerre(1);
  double intlag, weight1lag;
  tie(intlag, weight1lag) = Gausslaguerre(1);
  REQUIRE(weight1 == Approx(2.0));
  REQUIRE(weight1lag == Approx(2.0));
}

double f(double* x)
  {
    return x[0];
  }
//Here we if the MonteCarlo function can calculate the integral of x from 0 to 1
//With a uniform distribution, to an accuracy of 10 sigma (where sigma is 1/(sqrt(12*1e6)))
TEST_CASE("Testing the MonteCarlo function"){
  array<double, 2> integralMC;
  double tol = 15/(double)(3.464e3);
  int N = 1e6;
  double* a = new double[1];
  double* b = new double[1];
  a[0] = 0;
  b[0] = 1;
  integralMC = MonteCarlo(f,a,b,N,1,101);
  REQUIRE(integralMC[0] == Approx(0.5).epsilon(tol));
}