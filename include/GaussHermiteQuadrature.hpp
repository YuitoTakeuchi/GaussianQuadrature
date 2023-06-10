#pragma once
#include <cmath>

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>

template<int N>
class GaussHermiteQuadrature {
private:
    Eigen::VectorXd points = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);

public:
    // コンストラクタ
    GaussHermiteQuadrature () {
        // 行列を作る
        // pj = (aj x + bj)p_j-1 - cj p_j-2のとき
        // Hermiteでは aj = 1, bj = 0, cj = j-1
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N, N);
        double prev_a;//, prev_b, prev_c;
        prev_a = 1.0;
        // prev_b = 0.0;
        // prev_c = 0.0;
        for(int i = 1; i < N; ++i) {
            double a, b, c;
            a = 1.0;
            b = 0.0;
            c = 1.0*i;
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

    double integrate(double (*func)(double)) const {
        double ret = 0.0;
        for(int i = 0; i < N; ++i) {
            // 漸化式のところをすべて2倍すればsqrt(2)はなくなる
            // 理由は要検証
            ret += func(points(i) / sqrt(2)) * weights(i); 
        }
        return ret;
    }

    template<int Dimension>
    double integrage(double (*func)(Eigen::VectorXd x)) {
        double ret = 0.0;
        const long prod = pow(N, Dimension); // the number of iteration
        std::vector<int> product_set(Dimension, 0); // the cartesian product set of indices
        product_set[0] = -1;
        for(int i = 0; i < prod; ++i) {
            // create index list
            for(int idx = 0; idx < Dimension; ++idx) {
                if(++product_set[idx] < N) break;
                product_set[idx] = 0;
            }
            Eigen::VectorXd point = Eigen::VectorXd::Zero(Dimension);

            // calcurate i-th point and weight
            int cnt = -1; // index of point
            double weight = 1.0;
            for(auto &j: product_set) {
                point(++cnt) = points[j];
                weight *= weights[j];
            }
            ret += func(point / sqrt(2)) * weight;
        }
        return ret;
    }
};