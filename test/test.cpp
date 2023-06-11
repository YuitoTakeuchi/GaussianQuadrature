#include "GaussQuadrature.hpp"
#include <gtest/gtest.h>

namespace {
    class QuadratureTest: public :: testing::Test{}; 

    double test_func_1(double t) {
        return 1.0;
    }
    Eigen::VectorXd d_test_func_1(Eigen::VectorXd t) {
        Eigen::VectorXd ret = Eigen::VectorXd::Zero(t.rows());
        for(int i = 0; i < t.rows(); ++i) {
            ret(i) = 0.0;
        }
        return ret;
    }

    double test_func_2(double t) {
        return 3.0*t*t;
    }
    Eigen::VectorXd d_test_func_2(Eigen::VectorXd t) {
        Eigen::VectorXd ret = Eigen::VectorXd::Zero(t.rows());
        for(int i = 0; i < t.rows(); ++i) {
            ret(i) = 6.0*t(i);
        }
        return ret;
    }

    double test_func_3(double t) {
        return cos(t);
    }
    Eigen::VectorXd d_test_func_3(Eigen::VectorXd t) {
        Eigen::VectorXd ret = Eigen::VectorXd::Zero(t.rows());
        for(int i = 0; i < t.rows(); ++i) {
            ret(i) = -sin(t(i));
        }
        return ret;
    }

    Eigen::VectorXd vector_wrapper(double (*func)(double), Eigen::VectorXd t) {
        Eigen::VectorXd ret = Eigen::VectorXd::Zero(t.rows());
        for(int i = 0; i < t.rows(); ++i) {
            ret(i) = func(t(i));
        }
        return ret;
    }

    GaussLegendreQuadrature<1> glq1;
    TEST_F(QuadratureTest, integral_test_1_order_1) {EXPECT_NEAR(glq1.integrate(test_func_1), 2.0, 1e-12);}
    // TEST_F(QuadratureTest, integral_test_1_order_2) {EXPECT_NEAR(glq1.integrate(test_func_2), 2.0, 1e-12);}
    // TEST_F(QuadratureTest, integral_test_1_order_3) {EXPECT_NEAR(glq1.integrate(test_func_3), 2.0*sin(1.0), 1e-12);}
    GaussLegendreQuadrature<2> glq2;
    TEST_F(QuadratureTest, integral_test_2_order_1) {EXPECT_NEAR(glq2.integrate(test_func_1), 2.0, 1e-12);}
    TEST_F(QuadratureTest, integral_test_2_order_2) {EXPECT_NEAR(glq2.integrate(test_func_2), 2.0, 1e-12);}
    // TEST_F(QuadratureTest, integral_test_2_order_3) {EXPECT_NEAR(glq2.integrate(test_func_3), 2.0*sin(1.0), 1e-12);}
    GaussLegendreQuadrature<19> glq19;
    TEST_F(QuadratureTest, integral_test_19_order_1) {EXPECT_NEAR(glq19.integrate(test_func_1), 2.0, 1e-12);}
    TEST_F(QuadratureTest, integral_test_19_order_2) {EXPECT_NEAR(glq19.integrate(test_func_2), 2.0, 1e-12);}
    TEST_F(QuadratureTest, integral_test_19_order_3) {EXPECT_NEAR(glq19.integrate(test_func_3), 2.0*sin(1.0), 1e-12);}
    GaussLegendreQuadrature<20> glq20;
    TEST_F(QuadratureTest, integral_test_20_order_1) {EXPECT_NEAR(glq20.integrate(test_func_1), 2.0, 1e-12);}
    TEST_F(QuadratureTest, integral_test_20_order_2) {EXPECT_NEAR(glq20.integrate(test_func_2), 2.0, 1e-12);}
    TEST_F(QuadratureTest, integral_test_20_order_3) {EXPECT_NEAR(glq20.integrate(test_func_3), 2.0*sin(1.0), 1e-12);}

    class DerivativeTest: public::testing::Test{};
    TEST_F(DerivativeTest, integral_test_1_order_1) {
        Eigen::VectorXd expected = d_test_func_1(glq1.get_points());
        Eigen::VectorXd calculated = glq1.get_derivative(test_func_1);

        for(int i = 0; i < 1; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }

    TEST_F(DerivativeTest, derivative_test_2_order_1) {
        Eigen::VectorXd expected = d_test_func_1(glq2.get_points());
        Eigen::VectorXd calculated = glq2.get_derivative(test_func_1);

        for(int i = 0; i < 2; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }
    TEST_F(DerivativeTest, derivative_test_2_order_2) {
        Eigen::VectorXd expected = d_test_func_2(glq2.get_points());
        Eigen::VectorXd calculated = glq2.get_derivative(test_func_2);// * vector_wrapper(&test_func_2, glq2.get_points());

        for(int i = 0; i < 2; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }

    TEST_F(DerivativeTest, derivative_test_19_order_1) {
        Eigen::VectorXd expected = d_test_func_1(glq19.get_points());
        Eigen::VectorXd calculated = glq19.get_derivative(test_func_1);// * vector_wrapper(&test_func_1, glq19.get_points());

        for(int i = 0; i < 19; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }
    TEST_F(DerivativeTest, derivative_test_19_order_2) {
        Eigen::VectorXd expected = d_test_func_2(glq19.get_points());
        Eigen::VectorXd calculated = glq19.get_derivative(test_func_2);// * vector_wrapper(&test_func_2, glq19.get_points());

        for(int i = 0; i < 19; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }
    TEST_F(DerivativeTest, derivative_test_19_order_3) {
        Eigen::VectorXd expected = d_test_func_3(glq19.get_points());
        Eigen::VectorXd calculated = glq19.get_derivative(test_func_3);// * vector_wrapper(&test_func_3, glq19.get_points());

        for(int i = 0; i < 19; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }


    TEST_F(DerivativeTest, derivative_test_20_order_1) {
        Eigen::VectorXd expected = d_test_func_1(glq20.get_points());
        Eigen::VectorXd calculated = glq20.get_derivative(test_func_1);// * vector_wrapper(&test_func_1, glq20.get_points());

        for(int i = 0; i < 20; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }
    TEST_F(DerivativeTest, derivative_test_20_order_2) {
        Eigen::VectorXd expected = d_test_func_2(glq20.get_points());
        Eigen::VectorXd calculated = glq20.get_derivative(test_func_2);// * vector_wrapper(&test_func_2, glq20.get_points());

        for(int i = 0; i < 20; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }
    TEST_F(DerivativeTest, derivative_test_20_order_3) {
        Eigen::VectorXd expected = d_test_func_3(glq20.get_points());
        Eigen::VectorXd calculated = glq20.get_derivative(test_func_3);// * vector_wrapper(&test_func_3, glq20.get_points());

        for(int i = 0; i < 20; ++i) {
            EXPECT_NEAR(calculated(i), expected(i), 1e-12);
        }
    }


}