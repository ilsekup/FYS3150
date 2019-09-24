#include "catch.hpp"
#include "Project2.h"

TEST_CASE("Testing Eigenvalues"){
  int N = 4;
  bool potential = false;
  double max = 1.0;
  mat A = initialize(N, max, potential);
  mat R(N,N,fill::eye);
  iterative(A, R, N);
  vec eigen = get_eigvals(A, N);
  REQUIRE(eigen(0) ==Approx(6.1115));
  REQUIRE(eigen(1) == Approx(22.1115));
  REQUIRE(eigen(2) == Approx(41.8885));
  REQUIRE(eigen(3) == Approx(57.8885));
}

TEST_CASE("Testing Eigenvector Orthogonality"){
  int N = 4;
  bool potential = false;
  double max = 1.0;
  mat A = initialize(N, max, potential);
  mat R(N,N,fill::eye);
  iterative(A, R, N);
  //dot product of vector 1 and vector 2 should be zero
  double dot01 = R(0,0)*R(1,0) +R(0,1)*R(1,1) + R(0,2)*R(1,2) + R(0,3)*R(1,3);
  double dot12 = R(1,0)*R(2,0) +R(1,1)*R(2,1) + R(1,2)*R(2,2) + R(1,3)*R(2,3);
  double dot20 = R(2,0)*R(0,0) +R(2,1)*R(0,1) + R(2,2)*R(0,2) + R(2,3)*R(0,3);
  //dot product of vector 1 with vector 1 should be one
  double dot00 = R(0,0)*R(0,0) + R(0,1)*R(0,1) + R(0,2)*R(0,2) + R(0,3)*R(0,3);
  double dot22 = R(2,0)*R(2,0) + R(2,1)*R(2,1) + R(2,2)*R(2,2) + R(2,3)*R(2,3);
  REQUIRE(dot01 == Approx(0.000).epsilon(0.01));
  REQUIRE(dot12 == Approx(0.000).epsilon(0.01));
  REQUIRE(dot20 == Approx(0.000).epsilon(0.01));
  REQUIRE(dot00 == Approx(1.000).epsilon(0.01));
  REQUIRE(dot22 == Approx(1.000).epsilon(0.01));
}

