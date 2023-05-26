#include "GaussQuadrature.hpp"
#include <cmath>

double f(double x) {
    return x*x;
}

double f_o(double x) {
    return f(x) / sqrt(1-x*x);
}

int main() {
    int n;
    std::cin >> n;
    GaussQuadrature<LegendrePolynomial> gq;
    std::cout << gq.integrate(n, -0.5, 0.5, &f_o) << std::endl;
}