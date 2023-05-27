#include "GaussQuadrature.hpp"
#include <cmath>

double x(double t) {
    return cos(t);
}

Eigen::VectorXd x(Eigen::VectorXd t) {
    Eigen::VectorXd ret = t;
    for(int i = 0; i < ret.rows(); ++i) ret(i) = x(t(i));
    return ret;
}

int main() {
    const int n = 7;
    // GaussQuadrature<LegendrePolynomial> gq;
    // std::cout << gq.integrate(n, -0.5, 0.5, &f_o) << std::endl;
    // Eigen::VectorXd weights, points;

    GaussLegendreQuadrature<n> glq;
    std::cout << glq.derivative << std::endl << std::endl;
    std::cout << glq.points << std::endl << std::endl;

    Eigen::VectorXd D = glq.derivative * x(glq.points);

    std::cout << D << std::endl;
}