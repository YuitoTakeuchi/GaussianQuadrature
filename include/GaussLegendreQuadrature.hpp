#pragma once
#include <cmath>

#include <Eigen/Dense>
#include <lapacke.h>
#include <vector>
#include <iostream>
#include <iomanip>

// \int_-1^1 f(x)dx
// N is the number of nodes
template<int N>
class GaussLegendreQuadrature {
private:
    Eigen::VectorXd points = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd derivative_matrix = Eigen::MatrixXd::Zero(N, N+1);


public:
    // コンストラクタ
    GaussLegendreQuadrature () {
        // 行列を作る
        double d[N], e[N-1], z[N*N];
        memset(d, 0, N);

        double prev_a; //, prev_b, prev_c;
        prev_a = 1.0;
        // prev_b = 0.0;
        // prev_c = 0.0;
        for(int i = 1; i < N; ++i) {
            double a, b, c;
            a = (2.0*i+1.0)/(i+1.0);
            b = 0.0;
            c = i/(i+1.0);
            e[i-1] = sqrt(c/a/prev_a);
            prev_a = a;
            // prev_b = b;
            // prev_c = c;
        }

        LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', N, d, e, z, N);

        // 固有値と固有ベクトルから重みとノードを求める．
        for(int i = 0; i < N; ++i) {
            points(i)  = d[i];
            weights(i) = z[i*N] * z[i*N] * 2.0;
        }

        // 微分行列を計算
        // -1.0のところは適当なので（かぶらなければ）なんでもok
        for(int k = 0; k < N; ++k) {
            for(int l = 0; l < N; ++l) {
                if(k == l) {
                    for(int m = 0; m < N; ++m) {
                        if(m == k) continue;
                        derivative_matrix(k, l) += 1.0/(points(k) - points(m));
                    }
                    derivative_matrix(k, l) += 1.0 / (points(k) + 1.0);
                } else {
                    derivative_matrix(k, l) = 1.0/(points(l) - points(k));
                    for(int m = 0; m < N; ++m) {
                        if(m == l || m == k) continue;
                        derivative_matrix(k, l) *= (points(k) - points(m))/(points(l) - points(m));
                    }
                    derivative_matrix(k, l) *= (points(k) + 1.0)/(points(l) + 1.0);
                }
            }
            derivative_matrix(k, N) = 1.0/(-1.0 - points(k));
            for(int m = 0; m < N; ++m) {
                if(m == k) continue;
                derivative_matrix(k, N) *= (points(k) - points(m))/(-1.0 - points(m));
            }
        }
    }
    const Eigen::VectorXd& get_points() {return points;}
    const Eigen::VectorXd& get_weights() {return weights;}
    const Eigen::MatrixXd& get_derivative_matrix() {return derivative_matrix;}

    Eigen::VectorXd get_derivative(double (*func)(double)) const {
        Eigen::VectorXd _val = Eigen::VectorXd::Zero(N+1);
        for(int i = 0; i < N; ++i) _val(i) = func(points(i));
        _val(N) = func(-1.0);
        Eigen::VectorXd derivative = derivative_matrix * _val;
        return derivative;
    }

    double integrate(double (*func)(double)) {
        double ret = 0.0;
        for(int i = 0; i < N; ++i) {
            ret += func(points(i)) * weights(i);
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
            ret += func(point) * weight;
        }

        return ret;
    }
};