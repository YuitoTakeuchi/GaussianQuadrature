#include "GaussQuadrature.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

double x(double t) {
    return cos(t);
}

Eigen::VectorXd x(Eigen::VectorXd t) {
    const int n = t.rows();
    Eigen::VectorXd ret = Eigen::VectorXd::Zero(n+1);
    for(int i = 0; i < n; ++i) ret(i) = x(t(i));
    ret(n) = cos(1.0);
    return ret;
}

double f(double x) {
    double y =  2.0 + sqrt(2) * 0.2 * x;
    return cos(y*y);
}

int main() {
    constexpr int n = 20;

    // GaussLegendreQuadrature<n> glq;

    // Eigen::VectorXd d = std::move(glq.get_derivative(x));

    // std::cout << "calculated |expected   \n";
    // std::cout << "-----------|-----------\n";
    // for(int i = 0; i < n; ++i) {
    //     std::cout << std::setprecision(6) << std::setw(10);
    //     std::cout << d(i) << " | " << -sin(glq.get_points()(i)) << "\n";
    // }
    // std::cout << "-----------|-----------\n";

    GaussHermiteQuadrature <n> ghq;
    GaussHermiteQuadratureT<n> ghqt; // test for different coefficients


    // std::cout << "point      |weight     \n";
    // std::cout << "-----------|-----------\n";
    // for(int i = 0; i < n; ++i) {
    //     std::cout << std::setprecision(6) << std::setw(10);
    //     std::cout << ghq.get_points()(i) << " | " << ghq.get_weights()(i) << "\n";
    // }
    // std::cout << "-----------|-----------\n";

    double ret = ghq.integrate(f) / sqrt(M_PI);
    double retT = ghqt.integrate(f) / sqrt(M_PI);
    std::cout << "He (physits)         | H (probabilists)\n";
    std::cout << "---------------------|----------------------\n";
    std::cout << std::setprecision(15) << std::setw(20);
    std::cout << ret  << " | " << retT << std::endl;
    std::cout << "--------------------------------------------\n";
    std::cout << "ratio: " << ret / retT << std::endl;
    //     std::cout << std::setprecision(6) << std::setw(10);
    //     std::cout << ghq.get_points()(i) << " | " << ghq.get_weights()(i) << "\n";
}