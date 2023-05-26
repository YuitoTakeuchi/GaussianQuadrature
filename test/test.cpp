#include "GaussQuadrature.hpp"
#include <gtest/gtest.h>

namespace {
  class MyLibTest : public ::testing::Test{}; 

  TEST_F(MyLibTest, SqrFuncCheck) {
    EXPECT_EQ((3.0) * (3.0), 9.0);
    EXPECT_EQ((4.0) * (4.0), 16.0);
  }
}