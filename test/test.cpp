#include "GaussQuadrature.hpp"
#include <gtest/gtest.h>

namespace {
    class QuadratureTest: public :: testing::Test{}; 

    double test_func_1(double t) {
        return 1.0;
    }

    double test_func_2(double t) {
        return 3.0*t*t;
    }

    double test_func_3(double t) {
        return cos(t);
    }

    GaussLegendreQuadrature<1> glq1;
    TEST_F(QuadratureTest, integral_test_1_order_1) {EXPECT_DOUBLE_EQ(glq1.integrate(test_func_1), 2.0);}
    // TEST_F(QuadratureTest, integral_test_1_order_2) {EXPECT_DOUBLE_EQ(glq1.integrate(test_func_2), 2.0);}
    // TEST_F(QuadratureTest, integral_test_1_order_3) {EXPECT_DOUBLE_EQ(glq1.integrate(test_func_3), 2.0*sin(1.0));}
    GaussLegendreQuadrature<2> glq2;
    TEST_F(QuadratureTest, integral_test_2_order_1) {EXPECT_DOUBLE_EQ(glq2.integrate(test_func_1), 2.0);}
    TEST_F(QuadratureTest, integral_test_2_order_2) {EXPECT_DOUBLE_EQ(glq2.integrate(test_func_2), 2.0);}
    // TEST_F(QuadratureTest, integral_test_2_order_3) {EXPECT_DOUBLE_EQ(glq2.integrate(test_func_3), 2.0*sin(1.0));}
    GaussLegendreQuadrature<19> glq19;
    TEST_F(QuadratureTest, integral_test_19_order_1) {EXPECT_DOUBLE_EQ(glq19.integrate(test_func_1), 2.0);}
    TEST_F(QuadratureTest, integral_test_19_order_2) {EXPECT_DOUBLE_EQ(glq19.integrate(test_func_2), 2.0);}
    TEST_F(QuadratureTest, integral_test_19_order_3) {EXPECT_DOUBLE_EQ(glq19.integrate(test_func_3), 2.0*sin(1.0));}
    GaussLegendreQuadrature<20> glq20;
    TEST_F(QuadratureTest, integral_test_20_order_1) {EXPECT_DOUBLE_EQ(glq20.integrate(test_func_1), 2.0);}
    TEST_F(QuadratureTest, integral_test_20_order_2) {EXPECT_DOUBLE_EQ(glq20.integrate(test_func_2), 2.0);}
    TEST_F(QuadratureTest, integral_test_20_order_3) {EXPECT_DOUBLE_EQ(glq20.integrate(test_func_3), 2.0*sin(1.0));}
}