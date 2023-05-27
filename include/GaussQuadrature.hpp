#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

// ガウス求積
template<class OrthogonalPolynomial>
class GaussQuadrature {
private:
    OrthogonalPolynomial base_polynomial;
public:
    // n: 近似点の数が2n-1になる
    // a, b: 積分区間を[a, b]にとる
    // mu: 重み関数の積分． \int_a^b \omega(x) dx
    void calc_nodes(int n, Eigen::VectorXd& weights, Eigen::VectorXd& points, double a=-1, double b=1, double mu=2) {
        // 行列を作る
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
        double prev_a, prev_b, prev_c;
        base_polynomial.reccurence_coeffs(1, prev_a, prev_b, prev_c);
        for(int i = 1; i < n; ++i) {
            double a, b, c;
            base_polynomial.reccurence_coeffs(i+1, a, b, c);
            J(i-1, i-1) = -prev_b/prev_a;
            J(i-1, i) = sqrt(c/a/prev_a);
            J(i, i-1) = J(i-1, i);
            prev_a = a;
            prev_b = b;
            prev_c = c;
        }
        J(n-1, n-1) = -prev_b/prev_a;

        // 固有値と固有ベクトルを求める
        Eigen::EigenSolver<Eigen::MatrixXd> eigens(J);
        points = Eigen::VectorXd::Zero(n);
        weights = Eigen::VectorXd::Zero(n);
        Eigen::MatrixXcd eigen_vectors = eigens.eigenvectors();
        Eigen::VectorXcd eigen_values  = eigens.eigenvalues();

        // 固有値と固有ベクトルから重みとノードを求める．
        for(int i = 0; i < n; ++i) {
            Eigen::VectorXd eigen_vector =  eigen_vectors.col(i).real().normalized();
            points(i)  = 0.5*(b-a)*eigen_values(i).real() + (a+b)*0.5;
            weights(i) = mu*eigen_vector(0)*eigen_vector(0);
        }

        std::cout << points << std::endl;
    }

    double integrate(int n, double a, double b, double (*target_func)(double)) {
        double ans = 0;
        Eigen::VectorXd weights, points;
        calc_nodes(n, weights, points, a, b, base_polynomial.calc_mu(a, b));
        for(int i = 0; i < n; ++i) {
            ans += weights(i) * target_func(points(i));
        }
        return ans;
    }
};

class LegendrePolynomial {
public:
    // returns coefficients a_n, b_n, c_n
    // which satisfy p_n(x) = (a_nx + b_n)p_{n-1} (x) - c_n p_{n-2}(x)
    void reccurence_coeffs(int n, double& a, double& b, double& c) {
        a = (2.0*n-1.0)/n;
        b = 0.0;
        c = (n-1.0)/n;
    }

    // 重み関数の積分
    // \int_a^b w(x) dx
    // Legendre多項式の場合はw(x) = 1
    double calc_mu(double a, double b) {
        return b - a;
    }

};

class ChebyshevPolynomial {
    // 動作確認できていない
public:
    // returns coefficients a_n, b_n, c_n
    // which satisfy p_n(x) = (a_nx + b_n)p_{n-1} (x) - c_n p_{n-2}(x)
    void reccurence_coeffs(int n, double& a, double& b, double& c) {
        a = 2.0;
        b = 0.0;
        c = 1.0;
    }

    // 重み関数の積分
    // \int_a^b w(x) dx
    // Legendre多項式の場合はw(x) = 1
    double calc_mu(double a, double b) {
        return (asin(b) - asin(a) + b*sqrt(1-b*b) - a*sqrt(1-a*a))*0.5;
    }

};


// \int_-1^1 f(x)dx
// N is the number of nodes
template<int N>
class GaussLegendreQuadrature {
private:
public:
    Eigen::VectorXd points = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd derivative = Eigen::MatrixXd::Zero(N, N+1);
    GaussLegendreQuadrature () {
        // 行列を作る
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N, N);
        double prev_a, prev_b, prev_c;
        prev_a = 1.0;
        prev_b = 0.0;
        prev_c = 0.0;
        for(int i = 1; i < N; ++i) {
            double a, b, c;
            a = (2.0*i+1.0)/(i+1.0);
            b = 0.0;
            c = i/(i+1.0);
            J(i-1, i) = sqrt(c/a/prev_a);
            J(i, i-1) = J(i-1, i);
            prev_a = a;
            prev_b = b;
            prev_c = c;
        }

        // 固有値と固有ベクトルを求める
        Eigen::EigenSolver<Eigen::MatrixXd> eigens(J);
        points = Eigen::VectorXd::Zero(N+1);
        weights = Eigen::VectorXd::Zero(N+1);
        Eigen::MatrixXd eigen_vectors = eigens.eigenvectors().real();
        Eigen::VectorXd eigen_values  = eigens.eigenvalues().real();

        // 固有値と固有ベクトルから重みとノードを求める．
        std::vector<int> indices(N);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&eigen_values](int i1, int i2) {return eigen_values(i1) < eigen_values(i2);});
        int cnt = 1;
        points(0) = -1.0;
        for(int& i: indices) {
            points(cnt)  = eigen_values(i);
            Eigen::VectorXd eigen_vector =  eigen_vectors.col(i).normalized();
            weights(cnt) = 2.0*eigen_vector(0)*eigen_vector(0);
            cnt += 1;
        }

        for(int k = 1; k <= N; ++k) {
            for(int l = 0; l <= N; ++l) {
                if(k == l) {
                    for(int m = 0; m <= N; ++m) {
                        if(m == k) continue;
                        derivative(k-1, l) += 1.0/(points(k) - points(m));
                    }
                } else {
                    derivative(k-1, l) = 1.0/(points(l) - points(k));
                    for(int m = 0; m <= N; ++m) {
                        if(m == l || m == k) continue;
                        derivative(k-1, l) *= (points(k) - points(m))/(points(l) - points(m));
                    }
                }
            }
        }
    }
};

// // dynamic
// template<>
// class GaussLegendreQuadrature<-1> {

// }

// using GaussLegendreQuadrature = GuassLegendreQuadrature<-1>;