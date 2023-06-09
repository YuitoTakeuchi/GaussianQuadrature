#include "GaussQuadrature.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

double x(double t) {
    return cos(t);
}

Eigen::VectorXd x(Eigen::VectorXd t) {
    Eigen::VectorXd ret = t;
    for(int i = 0; i < ret.rows(); ++i) ret(i) = x(t(i));
    return ret;
}

double f(double x) {
    double y =  2.0 + sqrt(2) * 0.2 * x;
    return cos(y*y);
}

int main() {
    const int n = 30;
    // GaussQuadrature<LegendrePolynomial> gq;
    // std::cout << gq.integrate(n, -0.5, 0.5, &f_o) << std::endl;
    // Eigen::VectorXd weights, points;

    // GaussLegendreQuadrature<n> glq;
    // std::cout << glq.get_derivative() << std::endl << std::endl;
    // std::cout << glq.get_points() << std::endl << std::endl;

    // Eigen::VectorXd D = glq.get_derivative() * x(glq.get_points());

    // std::cout << D << std::endl;

    GaussHermiteQuadrature<n> ghq;


    std::cout << "point      |weight     \n";
    std::cout << "-----------|-----------\n";
    for(int i = 0; i < n; ++i) {
        std::cout << std::setprecision(6) << std::setw(10);
        std::cout << ghq.get_points()(i) << " | " << ghq.get_weights()(i) << "\n";
    }
    std::cout << "-----------|-----------\n";

    double ret = ghq.integrate(f) / sqrt(M_PI);
    std::cout << std::setprecision(15) << ret  << std::endl;
}