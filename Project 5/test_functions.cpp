#include "catch.hpp"
#include "PDEsolver.h"
#include "lib.h"

TEST_CASE("Testing the initial conditions implicit scheme"){
  int n = 10;
  int t_steps = 5;
  double dx = 1/(double) (n+1);
  double dt = 0.5*dx*dx;
  for(int i = 0; i < t_steps; i++){
    double *vin = implicit(n, dt, t_steps); //checking that the boundary conditions hold for all different timesteps
    REQUIRE(vin[0] == Approx(0));
    REQUIRE(vin[n+1] == Approx(1));
  }
}