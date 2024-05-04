// Only needed to eliminate intelliSennse givng false errors about Eigen.
// Not needed to compile, just gets rid of the false error squiggles
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include "Eigen/Dense"
#include <cmath>
#include <iostream>
#include <chrono>
#include "matrix_reader.h"
#define PI 3.1415926

using namespace std;
using namespace Eigen;

// Accepts any type of dense matrix. No sparse
template <typename Derived>
void SCludcmp(MatrixBase<Derived> &a, long n, int *indx) {
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *scaling;
    double TINY = 1.0e-20;

    scaling = new double[n];
    // *d = 1.0;
    for (i = 0; i < n; i++) {
        big = a.row(i).cwiseAbs().maxCoeff(); // the max abs(value) of i_th row

        if (big == 0.0) {
            printf("Singular matrix in routine ludcmp");
            exit(1);
        }
        scaling[i] = 1.0 / big;
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            sum = a(i, j);
            for (k = 0; k < i; k++)
                sum -= a(i, k) * a(k, j);

            a(i, j) = sum;
        }
        big = 0.0;
        for (i = j; i < n; i++) {
            sum = a(i, j);
            for (k = 0; k < j; k++)
                sum -= a(i, k) * a(k, j);

            a(i, j) = sum;
            // if(j == n-1) // print a(n,n)
            cout << "a_ij " << i << "," << j << " " << a(i, j) << endl;

            if ((dum = scaling[i] * abs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        // cout << "imax " << imax << endl;
        if (j != imax) {
            a.row(imax).swap(a.row(j));
            scaling[imax] = scaling[j];
        }
        indx[j] = imax;
        if (a(j, j) == 0.0)
            a(j, j) = TINY;
        if (j != n - 1) {
            dum = (1.0 / a(j, j));
            for (i = j + 1; i < n; i++)
                a(i, j) *= dum;
        }
    }
}

template <typename Derived>
void backsub(MatrixBase<Derived> &a, int n, int *indx, double b[]) {
    int i, ii = 0, ip, j;
    double sum;
    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii != 0)
            for (j = ii - 1; j < i; j++)
                sum -= a(i, j) * b[j];
        else if (sum != 0.0)
            ii = i + 1;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < n; j++)
            sum -= a(i, j) * b[j];
        b[i] = sum / a(i, i);
    }
}

int main(int argc, char* argv[]) {
    // double b[] = {2 * PI, 5 * PI, -8 * PI};
    //double b[] = {1, 1, -2, -2};

     MatrixXf mat;
    // load the matrix from the arguments
    if (argc != 2) {
        cout << "Usage: ./PCludcmp filename.mtx NUM_THREADS" << endl;
        cout << "CSV or matrix market .mtx files are allowed, file extension "
                "must be correct"
             << endl;
        return 1;
    }
    string filename = argv[1];
    //cout << "Threads: " << stoi(argv[2]) << endl;  
    string ext = filename.substr(filename.find_last_of(".") + 1);
    transform(ext.begin(), ext.end(), ext.begin(),
              [](unsigned char c) { return std::tolower(c); });
    if (ext == "csv") {
        mat = load_csv(filename);
    } else if (ext == "mtx") {
        mat = read_matrix_market(filename);
    } else {
        cout << "Usage: ./PCludcmp filename.mtx" << endl;
        cout << "CSV or matrix market .mtx files are allowed, file extension "
                "must be correct"
             << endl;
        cout << "Bad extension" << endl;
        return 1;
    }
    cout << "Starting timer" << endl;
    long n = mat.rows();
    int *indx = new int[n];
    // boost::timer myTimer;
    auto start = chrono::high_resolution_clock::now();
   // SCludcmp(mat, n, indx);
    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    cout << "Stopping timer" << endl;
    cout << "Time: " << duration.count() << endl;
    cout << mat.lu().matrixLU() << endl;
    //cout << mat << endl;
    return 0;
    // init matrix
    // a << 1.0, 2.0, -1.0, 6.0, -5.0, 4.0, -9.0, 8.0, -7.0;
    //a << 1, 0, 0, 1, 3, 2, 1, 0, 0, 0, 0, -1, 0, -1, 0, -2;
    // a = MatrixXd::Random(n, n);
    // // a = a + (MatrixXd::Constant(n, n, 1) *
    // //          100.0); // random n x n matrix doubles between -100, 100
    // //  cout << a << endl;
    // //  cout << a.lu().matrixLU() << endl; // print eigen's lu decomp partial
    // //  pivot
    // printf("\n");
    // // in place LU decomp from Eigen for testing.
    // // PartialPivLU<Ref<RowMajorMatrixXd>> lu(a);
    // // cout << a.col(1).maxCoeff(&index) << " index" << index << endl;
    // SCludcmp(a, n, indx);
    // cout << a << endl; // print our version
    // //backsub(a, n, indx, b);
    // cout << b << endl;

    return 0;
}