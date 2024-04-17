
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

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "To use this function, provide a <filename.csv> and <matrix_size>" << endl;
        cerr << "Usage: " << argv[0] << " <filename.csv>" << " <matrix_size>" << endl;
        return 1;
    }

    string filename = argv[1];
    int matrix_size = stoi(argv[2]);
    int** result = readCSV(filename, matrix_size);
    if (PRINT_MATRIX_BEFORE) {
        print_matrix(result, matrix_size);
    }

    free_dynamic_array(result, matrix_size);
    return 0;
}
