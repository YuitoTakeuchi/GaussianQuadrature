#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

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
        points = Eigen::VectorXd::Zero(N+1);
        weights = Eigen::VectorXd::Zero(N+1);
        Eigen::MatrixXd eigen_vectors = eigens.eigenvectors().real();
        Eigen::VectorXd eigen_values  = eigens.eigenvalues().real();

        // 固有値と固有ベクトルから重みとノードを求める．
        // 点を小さい順に並び替える
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

    Eigen::VectorXd get_points {return points;}
    Eigen::VectorXd get_weights {return weights;}
    Eigen::MatrixXd get_derivative {return derivative;}
};