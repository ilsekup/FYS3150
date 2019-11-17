#include "catch.hpp"
#include "isingmodel.h"
#include "lib.h"
//Testing the initialize function, if ordered is true all matrix elements should be 1
TEST_CASE("Testing the initialize function"){
  long seed = time(NULL);
  int **spin_matrix = (int**) matrix(2, 2, sizeof(int));
  double E = 0;
  double M = 0;
  initialize(2, 1.0, spin_matrix, E, M, seed, true);
  int oneone = spin_matrix[0][0];
  int onetwo = spin_matrix[0][1];
  int twoone = spin_matrix[1][0];
  int twotwo = spin_matrix[1][1];
  initialize(2, 1.0, spin_matrix, E, M, seed, false);
  REQUIRE(oneone == Approx(1));
  REQUIRE(onetwo == Approx(1));
  REQUIRE(twoone == Approx(1));
  REQUIRE(twotwo == Approx(1));
}