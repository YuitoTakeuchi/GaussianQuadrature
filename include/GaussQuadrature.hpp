#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <cmath>

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