#include "catch.hpp"
#include "Project2.h"
/*
TEST_CASE("Testing Maximum Values"){
  int N = 3;
  mat A = initialize(N);
  int *k;
  int *l;
  double max_A = maxoff(A, N, k, l);
  REQUIRE(*k==2);
  REQUIRE(*l==1);
}
*/
TEST_CASE("Testing Eigenvalues"){
  int N = 4;
  potential = false;
  max = 1;
  mat A = initialize(N, max, potential);
  mat R(N,N,fill::eye);
  iterative(A, R, N);
  vec eigen = get_eigvals(A, N);
  REQUIRE(eigen(0) ==Approx(6.1115));
  REQUIRE(eigen(1) == Approx(22.1115));
  REQUIRE(eigen(2) == Approx(41.8885));
  REQUIRE(eigen(3) == Approx(57.8885));
}


