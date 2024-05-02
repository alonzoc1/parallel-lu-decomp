#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include "Eigen/Dense"
#include <cmath>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _MPI
#include <mpi.h>
#endif
#define PI 3.1415926

using namespace std;
using namespace Eigen;

template <typename Derived>
void PCludcmp(MatrixBase<Derived> &a, long n, int *indx) {
    int i, imax, j, k;
    double big, dum, sum;
    double *scaling;
    double TINY = 1.0e-20;
    // mpi parameters
    int comm_sz, my_rank;
    // MPI_Init()

    scaling = new double[n];

    // implicit scaling. We find the scaling factor and save it.
    // no actual scaling is made here. Row i is implicitly scaled by scaling[i]
#pragma omp parallel for default(none) private(i, big) shared(scaling, a, n)
    for (i = 0; i < n; i++) {
        big = a.row(i).cwiseAbs().maxCoeff(); // the max abs(value) of i_th row
        if (big == 0.0) {
            printf("Singular matrix in routine ludcmp");
            exit(1);
        }
        scaling[i] = 1.0 / big; // save the scaling factor
    }

    // main crout's loop. iterate by column
    for (j = 0; j < n; j++) {
// build upper trianfular part up to j-1 row since eache thread will have its own set of i's,
// no overlap the a(k,j) term would have already been build in previous iterations of j
#pragma omp parallel for default(none) private(i, k, sum) shared(a, j)         \
    schedule(guided)
        for (i = 0; i < j; i++) {
            sum = a(i, j);
            for (k = 0; k < i; k++)
               
                sum -= a(i, k) * a(k, j); // we are not reducing to sum

            a(i, j) = sum;
        }
    
#pragma omp parallel for default(none) private(i, sum, k)                      \
    shared(a, j, n, big, imax, cout) schedule(static)
        for (i = j; i < n; i++) {
            sum = a(i, j);
            for (k = 0; k < j; k++)
                sum -= a(i, k) * a(k, j);

            a(i, j) = sum;
            // if(j == n-1) // print a(n,n)
            cout << "a_ij " << i << "," << j << " " << a(i, j) << endl;
        }

        // look for pivot. This area is not worth parallelizing since we will have to make the
        // conditional a critical section. We interested in the row-index of max(big), not it's value.
        // the operatio is scalar-scalar multiplication
        big = 0.0;
        for (i = j; i < n; i++) {
            dum = scaling[i] * abs(a(i, j));
            if (dum >= big) {
                big = dum;
                imax = i;
            }
        }

        // swap pivot rows. No parallelization
        if (j != imax) {
            a.row(imax).swap(a.row(j));
            scaling[imax] = scaling[j];
        }

        indx[j] = imax; // needed for backsubtitution if solveing Ax = b
        if (a(j, j) == 0.0)
            a(j, j) = TINY;
        if (j != n - 1) {
            dum = (1.0 / a(j, j));
#pragma omp parallel for default(none) private(i, k) shared(j, n, a, dum)      \
    schedule(static)
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

int main(int argc, char **argv) {
    long n = (long)stoi(argv[1]);
    MatrixXd a(n, n); // test matrix
    // initParallel();
    int *indx = new int[n];
    double b[] = {1, 1, -2, -2};
    // init matrix
    // a << 1.0, 2.0, -1.0, 6.0, -5.0, 4.0, -9.0, 8.0, -7.0;
    a << 1, 0, 0, 1, 3, 2, 1, 0, 0, 0, 0, -1, 0, -1, 0, -2;
    // a = MatrixXd::Random(n, n);
    // a = a + (MatrixXd::Constant(n, n, 1) *
    //          100.0); // random n x n matrix doubles between -100, 100
    //  cout << a << endl;
    //  cout << a.lu().matrixLU() << endl; // print eigen's lu decomp partial
    //  pivot
    printf("\n");
    // in place LU decomp from Eigen for testing.
    // PartialPivLU<Ref<RowMajorMatrixXd>> lu(a);
    // cout << a.col(1).maxCoeff(&index) << " index" << index << endl;
    PCludcmp(a, n, indx);
    std::cout << a << endl; // print our version
    backsub(a, n, indx, b);

    return 0;
}
