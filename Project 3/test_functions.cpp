#include "catch.hpp"
#include "Gauss.h"

TEST_CASE("Testing Eigenvalues"){
  double integralgauss = Gausslaguerre(10);
  REQUIRE(integralgauss ==Approx(0.19277).epsilon(0.1));
}