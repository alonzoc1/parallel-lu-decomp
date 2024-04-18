
#include <iostream>
#include <string>
#include <mpi.h>
#include "matrix_reader.h"

using namespace std;

// set to true to verify the matrix is being read in as expected
const bool PRINT_MATRIX_BEFORE = true;

void free_dynamic_array(int** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        delete[] arr[i];
    }
    delete[] arr;
}

void print_matrix(int** arr, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << arr[i][j] << "\t";
        }
        cout << endl;
    }
}

void lu_decomposition(int** data, int size) {
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "To use this function, provide a <filename.mtx>" << endl;
        cerr << "Usage: " << argv[0] << " <filename.mtx>" << endl;
        return 1;
    }
    string filename = argv[1];
    Eigen::MatrixXf mat = read_matrix_market(filename);
    cout << mat(1, 1) << endl;
    return 0;
}
