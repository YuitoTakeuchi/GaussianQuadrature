#include "GaussQuadrature.hpp"
#include <cmath>

double x(double t) {
    return 1;
}

Eigen::VectorXd x(Eigen::VectorXd t) {
    Eigen::VectorXd ret = t;
    for(int i = 0; i < ret.rows(); ++i) ret(i) = x(t(i));
    return ret;
}

int main() {
    const int n = 4;
    // GaussQuadrature<LegendrePolynomial> gq;
    // std::cout << gq.integrate(n, -0.5, 0.5, &f_o) << std::endl;
    // Eigen::VectorXd weights, points;

    GaussLegendreQuadrature<n> glq;
    std::cout << glq.derivative << std::endl << std::endl;
    std::cout << glq.points << std::endl << std::endl;

    Eigen::VectorXd D = glq.derivative * x(glq.points);


    // for(int i = 1; i <= n; ++i){
    //     double t = glq.points(i);
    //     double expected = 1;
    //     double a = 0.0;
    //     for(int j = 0; j <= n; ++j) {
    //         a += glq.derivative(i-1, j) * x(glq.points(j));
    //     }
    //     std::cout << "x: " << glq.points(i) << " a: " << a << " expected: " << expected << "\n";
    // }

    std::cout << D << std::endl;
}