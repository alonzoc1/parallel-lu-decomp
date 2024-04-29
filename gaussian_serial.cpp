
#include <iostream>
#include <string>
#include "matrix_reader.h"

using namespace std;
using namespace Eigen;

void lu_gaussian(MatrixXf mat, MatrixXf& L, MatrixXf& U) {
    unsigned long n = mat.rows();
    L = Eigen::MatrixXf::Identity(n, n);
    U = MatrixXf::Zero(n, n);
    for (int k = 0; k < n; ++k) {
        U(k, k) = mat(k, k);
        for (int i = k + 1; i < n; ++i) {
            L(i, k) = mat(i, k) / U(k, k);
            U(k, i) = mat(k, i);
        }
        for (int i = k + 1; i < n; ++i) {
            for (int j = k + 1; j < n; ++j) {
                mat(i, j) -= L(i, k) * U(k, j);
            }
        }
    }
}

int main(int argc, char* argv[]) {
    //Eigen::MatrixXf mat = read_matrix_market("data/mcca.mtx");
    MatrixXf mat = load_csv("data/mini.csv");
    long n = mat.rows();
    cout << "Input matrix:" << endl;
    cout << mat << endl;
    MatrixXf L, U;
    lu_gaussian(mat, L, U);
    cout << "L:" << endl;
    cout << L << endl;
    cout << "U:" << endl;
    cout << U << endl;
    cout << "L * U equals:" << endl;
    cout << (L*U) << endl;
    return 0;
}
