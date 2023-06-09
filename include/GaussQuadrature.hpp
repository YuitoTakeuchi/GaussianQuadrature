#pragma once
#include <cmath>

#include <Eigen/Dense>

// \int_-1^1 f(x)dx
// N is the number of nodes
template<int N>
class GaussLegendreQuadrature {
private:
    Eigen::VectorXd points = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd derivative = Eigen::MatrixXd::Zero(N, N+1);


public:
    // コンストラクタ
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
        points = Eigen::VectorXd::Zero(N);
        weights = Eigen::VectorXd::Zero(N);
        Eigen::MatrixXd eigen_vectors = eigens.eigenvectors().real();
        Eigen::VectorXd eigen_values  = eigens.eigenvalues().real();

        // 固有値と固有ベクトルから重みとノードを求める．
        for(int i = 0; i < N; ++i) {
            points(i)  = eigen_values(i);
            Eigen::VectorXd eigen_vector =  eigen_vectors.col(i).normalized();
            weights(i) = eigen_vector(0)*eigen_vector(0) * 2.0;
        }

        for(int k = 0; k < N; ++k) {
            for(int l = 0; l < N; ++l) {
                if(k == l) {
                    for(int m = 0; m < N; ++m) {
                        if(m == k) continue;
                        derivative(k, l) += 1.0/(points(k) - points(m));
                    }
                    derivative(k, l) += 1.0 / (points(k) + 1.0);
                } else {
                    derivative(k, l) = 1.0/(points(l) - points(k));
                    for(int m = 0; m < N; ++m) {
                        if(m == l || m == k) continue;
                        derivative(k, l) *= (points(k) - points(m))/(points(l) - points(m));
                    }
                    derivative(k, l) *= (points(k) + 1.0)/(points(l) + 1.0);
                }
            }
            // l == N
            derivative(k, N) = 1.0/(-1 - points(k));
            for(int m = 0; m < N; ++m) {
                if(m == k) continue;
                derivative(k, N) *= (points(k) - points(m))/(-1.0 - points(m));
            }
        }
    }
    Eigen::VectorXd& get_points() {return points;}
    Eigen::VectorXd& get_weights() {return weights;}
    Eigen::MatrixXd& get_derivative() {return derivative;}

    double integrate(double (*func)(double)) {
        double ret = 0.0;
        for(int i = 0; i < N; ++i) {
            ret += func(points(i)) * weights(i);
        }
        return ret;
    }
};

template<int N>
class GaussHermiteQuadrature {
private:
    Eigen::VectorXd points = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd derivative = Eigen::MatrixXd::Zero(N, N+1);


public:
    // コンストラクタ
    GaussHermiteQuadrature () {
        // 行列を作る
        // pj = (aj x + bj)p_j-1 - cj p_j-2のとき
        // Hermiteでは aj = 1, bj = 0, cj = j-1
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N, N);
        double prev_a;//, prev_b, prev_c;
        prev_a = 2.0;
        // prev_b = 0.0;
        // prev_c = 0.0;
        for(int i = 1; i < N; ++i) {
            double a, b, c;
            a = 2.0;
            b = 0.0;
            c = 2.0*i;
            J(i-1, i) = sqrt(c/a/prev_a);
            J(i, i-1) = J(i-1, i);
            prev_a = a;
            // prev_b = b;
            // prev_c = c;
        }

        // 固有値と固有ベクトルを求める
        Eigen::EigenSolver<Eigen::MatrixXd> eigens(J);
        points = Eigen::VectorXd::Zero(N);
        weights = Eigen::VectorXd::Zero(N);
        Eigen::MatrixXd eigen_vectors = eigens.eigenvectors().real();
        Eigen::VectorXd eigen_values  = eigens.eigenvalues().real();

        // 固有値と固有ベクトルから重みとノードを求める．
        for(int i = 0; i < N; ++i) {
            points(i)  = eigen_values(i);
            Eigen::VectorXd eigen_vector =  eigen_vectors.col(i).normalized();
            weights(i) = eigen_vector(0)*eigen_vector(0) * sqrt(M_PI);
        }
    }
    const Eigen::VectorXd& get_points() const {return points;}
    const Eigen::VectorXd& get_weights() const {return weights;}
    const Eigen::MatrixXd& get_derivative() const {return derivative;}

    double integrate(double (*func)(double)) const {
        double ret = 0.0;
        for(int i = 0; i < N; ++i) {
            ret += func(points(i)) * weights(i);
        }
        return ret;
    }
};