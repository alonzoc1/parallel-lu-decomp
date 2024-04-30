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

using namespace std;
using namespace Eigen;

template <typename Derived> void PCludcmp(MatrixBase<Derived> &a, long n) {
    int i, imax, j, k;
    double big, dum, sum;
    double *scaling;
    double TINY = 1.0e-20;

    scaling = new double[n];

    // implicit scaling. We find the scaling factor and save it.
    // no actual scaling is made here. Row i is implicitly scaled by scaling[i]
    // TODO: This part is easily doable with openMP
#pragma omp parallel for default(none) private(i, big) shared(scaling, a, n)
    for (i = 0; i < n; i++) {
        // big = 0.0;
        big = a.row(i).cwiseAbs().maxCoeff(); // the max abs(value) of i_th row
        // for (j = 0; j < n; j++){
        //     if(abs(a(i,j)) > big)
        //         big = abs(a(i,j));
        // }
        if (big == 0.0) {
            printf("Singular matrix in routine ludcmp");
            exit(1);
        }
        scaling[i] = 1.0 / big; // save the scaling factor
    }

    // main crout's loop. iterate by column
    for (j = 0; j < n; j++) {
        // build upper trianfular part up to j-1 row
#pragma omp parallel for default(none) private(i, k, sum) shared(a, j)
        for (i = 0; i < j; i++) { // this can be parallelize with Openmp
            sum = a(i, j);
            for (k = 0; k < i; k++)
                // since eache thread will have its own set of i's, no overlap
                // the a(k,j) term would have already been build in previous
                // iterations of j
                sum -= a(i, k) * a(k, j); // we are not reducing to sum, sum is
                                          // a private variable

            // no need for critical section because i's do not overlap
            a(i, j) = sum;
        }

        // find larget pivot element. Also compute i = j up to n rows
        big = 0.0;
        // parallelization here is very similar to the previous loop
        // FIXME: Problem seems to be with his loop. when I un-parallelize it it works
#pragma omp parallel for default(none) private(i, sum, k, dum)                 \
        shared(a, scaling, imax, j, n) reduction(max : big)
        for (i = j; i < n; i++) {
            sum = a(i, j);
            for (k = 0; k < j; k++)
                sum -= a(i, k) * a(k, j);

            a(i, j) = sum;
            // check if current pivot is the biggest. Only one thread at a time
            // prevents read-write racing condition
#pragma omp critical
            {
                if ((dum = scaling[i] * abs(sum)) >= big) {
                    big = dum; // return only the max value of big
                    imax = i;  // we are keeping track of the index of the pivot
                }
            }
        }

        // swap pivot rows.
        if (j != imax) {
            a.row(imax).swap(a.row(j));
            // for (k = 0; k < n; k++){
            // 	dum = a(imax,k);
            // 	a(imax,k) = a(j,k);
            // 	a(j,k) = dum;
            // }
            scaling[imax] = scaling[j];
        }

        if (a(j, j) == 0.0)
            a(j, j) = TINY;
        if (j != n - 1) {
            dum = (1.0 / a(j, j));
# pragma omp parallel for default(none) private(i) shared(j, n, a, dum)
            for (i = j + 1; i < n; i++)
                a(i, j) *= dum;
        }
    }
}

int main() {
    long n = 3;
    MatrixXd a(n, n); // test matrix

    // init matrix
    a << 1.0, 2.0, -1.0, 6.0, -5.0, 4.0, -9.0, 8.0, -7.0;
    cout << a << endl;

    // in place LU decomp from Eigen for testing.
    // PartialPivLU<Ref<RowMajorMatrixXd>> lu(a);

    PCludcmp(a, n);
    cout << a << endl;

    return 0;
}
