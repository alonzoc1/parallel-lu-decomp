#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;

void SCludcmp(MatrixXd& a, long n, long* indx) {
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;
    double TINY = 1.0e-20;

	vv = new double[n];
	// *d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++){
            temp = abs(a(i,j));
            if(temp > big) big = temp;
        }
		if (big == 0.0){ 
            printf("Singular matrix in routine ludcmp");
            exit(1);
        }
	    vv[i] = 1.0 / big;
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
			if ((dum = vv[i] * abs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
            // a.row(imax).swap(a.row(j));
			for (k = 0; k < n; k++){
				dum = a(imax,k);
				a(imax,k) = a(j,k);
				a(j,k) = dum;
			}
			// *d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
        if(a(j,j) == 0.0) a(j,j) = TINY;
		if (j != n-1) {
			dum = (1.0 / a(j,j));
			for (i = j + 1; i < n; i++) a(i,j) *= dum;
		}
	}
}

int main(){
    MatrixXd a(3,3);
    long* indx = new long[3];
    a << 1.0, 2.0, -1.0,
        6.0, -5.0, 4.0,
        -9.0, 8.0, -7.0;
        cout << a <<endl;

        SCludcmp(a, 3, indx);
        cout << a << endl;
    
    return 0;
}