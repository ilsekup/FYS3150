#include "catch.hpp"
#include "Gauss.h"

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