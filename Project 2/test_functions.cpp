#include "catch.hpp"
#include "Project2.h"

TEST_CASE("Testing Maximum Values"){
  int N = 3;
  mat A = initialize(N);
  int *k;
  int *l;
  double max_A = maxoff(A, N, k, l);
  REQUIRE(*k==2);
  REQUIRE(*l==1);
}
TEST_CASE("Testing Eigenvalues"){
  int N = 3;
  mat A = initialize(N);
  iterative(A, N);
  vec eigen = get_eigvals(A, N);
  REQUIRE(eigen(0) ==Approx(5.2721));
  REQUIRE(eigen(1) == Approx(18.0000));
  REQUIRE(eigen(2) == Approx(30.7279));
}


