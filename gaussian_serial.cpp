
#include <iostream>
#include <string>
#ifndef BOOST_TIMER_ENABLE_DEPRECATED
#define BOOST_TIMER_ENABLE_DEPRECATED
#endif
#include <boost/timer.hpp>
#include "matrix_reader.h"

using namespace std;
using namespace Eigen;

void lu_gaussian(MatrixXf mat, MatrixXf& L, MatrixXf& U) {
    unsigned long n = mat.rows();
    L = MatrixXf::Identity(n, n);
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
    MatrixXf mat;
    // load the matrix from the arguments
    if (argc != 2) {
        cout << "Usage: ./gaussian_serial filename.mtx" << endl;
        cout << "CSV or matrix market .mtx files are allowed, file extension must be correct" << endl;
        return 1;
    }
    string filename = argv[1];
    string ext = filename.substr(filename.find_last_of(".") + 1);
    transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return std::tolower(c); });
    if (ext == "csv") {
        mat = load_csv(filename);
    } else if (ext == "mtx") {
        mat = read_matrix_market(filename);
    } else {
        cout << "Usage: ./gaussian_serial filename.mtx" << endl;
        cout << "CSV or matrix market .mtx files are allowed, file extension must be correct" << endl;
        cout << "Bad extension" << endl;
        return 1;
    }
    cout << "Starting timer" << endl;
    boost::timer myTimer;
    long n = mat.rows();
    /*
    cout << "Input matrix:" << endl;
    cout << mat << endl;
    */
    MatrixXf L, U;
    lu_gaussian(mat, L, U);
    /*
    cout << "L:" << endl;
    cout << L << endl;
    cout << "U:" << endl;
    cout << U << endl;
    cout << "L * U equals:" << endl;
    cout << (L*U) << endl;
    */
    double timer = myTimer.elapsed();
    cout << "Stopping timer" << endl;
    MatrixXf result = L * U;
    MatrixXf reload = ext == "csv" ? load_csv(filename) : read_matrix_market(filename);
    cout << "Correct: " << ((result.isApprox(reload)) ? "Yes" : "No") << endl;
    cout << "Time: " << timer << endl;
    return 0;
}
