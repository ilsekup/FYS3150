#include "catch.hpp"
#include "PDEsolver.h"
#include "lib.h"
//Testing the initialize function, if ordered is true all matrix elements should be 1
TEST_CASE("Testing the initial conditions implicit scheme"){
  int n = 10;
  int t_steps = 5;
  double *vin = implicit(n, t_steps);
  REQUIRE(vin[0] == Approx(0));
  REQUIRE(vin[n] == Approx(1));
}