// Only needed to eliminate intelliSennse givng false errors about Eigen.
// Not needed to compile, just gets rid of the false error squiggles 
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <cmath>
#include <iostream>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

// Accepts any type of dense matrix. No sparse
template <typename Derived>
void SCludcmp(MatrixBase<Derived>& a, long n) {
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *scaling;
    double TINY = 1.0e-20;

	scaling = new double[n];
	// *d = 1.0;
	for (i = 0; i < n; i++) {
		big = a.row(i).cwiseAbs().maxCoeff(); // the max abs(value) of i_th row
		// big = 0.0;
		// for (j = 0; j < n; j++){
        //     temp = abs(a(i,j));
        //     if(temp > big) big = temp;
        // }
		if (big == 0.0){ 
            printf("Singular matrix in routine ludcmp");
            exit(1);
        }
	    scaling[i] = 1.0 / big;
    }
	
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a(i,j);
			for (k = 0; k < i; k++) 
                sum -= a(i,k) * a(k,j);
			
            a(i,j) = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a(i,j);
			for (k = 0; k < j; k++) 
                sum -= a(i,k) * a(k,j);
			
            a(i,j) = sum;
			if ((dum = scaling[i] * abs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
            a.row(imax).swap(a.row(j));
			// for (k = 0; k < n; k++){
			// 	dum = a(imax,k);
			// 	a(imax,k) = a(j,k);
			// 	a(j,k) = dum;
			// }
			scaling[imax] = scaling[j];
		}
        if(a(j,j) == 0.0) a(j,j) = TINY;
		if (j != n-1) {
			dum = (1.0 / a(j,j));
			for (i = j + 1; i < n; i++) a(i,j) *= dum;
		}
	}
}

int main(){
	long n = 3;
	MatrixXd a(n,n); // test matrix
	

	// init matrix
    a << 1.0, 2.0, -1.0,
        6.0, -5.0, 4.0,
        -9.0, 8.0, -7.0;
        cout << a <<endl;

		// in place LU decomp from Eigen for testing.
		//PartialPivLU<Ref<RowMajorMatrixXd>> lu(a);


        SCludcmp(a, n);
        cout << a << endl;


    
    return 0;
}