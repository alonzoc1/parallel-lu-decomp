
#include <iostream>
#include <string>
#include <mpi.h>
#include "matrix_reader.h"

using namespace std;
using namespace Eigen;

void lu_parallel(int rank, int procs, MatrixXf mat, int n) {
    for (int k = 0; k < n; ++k) {
        // broadcast current pivot row if its my turn
        MPI_Bcast(mat.row(k).data(), n, MPI_FLOAT, k % procs, MPI_COMM_WORLD);

        for (int i = k + 1; i < n; ++i) {
            if (i % procs == rank) {
                double factor = mat(i, k) / mat(k, k);
                for (int j = k; j < n; ++j) {
                    mat(i, j) -= factor * mat(k, j);
                }
                mat(i, k) = factor;
            }
        }
    }
}

void split_lu(MatrixXf& lu, MatrixXf& l, MatrixXf& u) {
    unsigned long n = lu.rows();
    l = MatrixXf::Zero(n, n);
    u = MatrixXf::Zero(n, n);
    // split
    for (unsigned long i = 0; i < n; ++i) {
        for (unsigned long j = 0; j < n; ++j) {
            if (i > j) {
                l(i, j) = lu(i, j); // lower
            } else {
                u(i, j) = lu(i, j); // upper
                if (i == j)
                    l(i, j) = 1.0; // l diagonal
            }
        }
    }
}

int main(int argc, char *argv[]) {
    // set up mpi
    MPI_Init(&argc, &argv);
    int rank, procs;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MatrixXf mat = MatrixXf::Ones(stol(argv[2]), stol(argv[2])); // initializes the memory space in each process
    unsigned long n = stoul(argv[2]);
    double timing;
    string ext;
    string filename;

    if (rank == 0) {
        // load the matrix from the arguments
        if (argc != 3) {
            cout << "Usage: ./gaussian_mpi filename.mtx n" << endl;
            cout << "CSV or matrix market .mtx files are allowed, file extension must be correct" << endl;
            cout << "n is the square size of the matrix (i.e. a 3x3 matrix is n=3)" << endl;
            return 1;
        }
        filename = argv[1];
        ext = filename.substr(filename.find_last_of(".") + 1);
        transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return std::tolower(c); });
        if (ext == "csv") {
            mat = load_csv(filename);
        } else if (ext == "mtx") {
            mat = read_matrix_market(filename);
        } else {
            cout << "Usage: ./gaussian_mpi filename.mtx" << endl;
            cout << "CSV or matrix market .mtx files are allowed, file extension must be correct" << endl;
            cout << "Bad extension" << endl;
            return 1;
        }
        // print input
        //cout << "Input mat: " << endl;
        //cout << mat << endl;
        cout << "Starting timer..." << endl;
        timing = MPI_Wtime();
    }
    // give data to processes
    MPI_Bcast(mat.data(), n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // run it!
    lu_parallel(rank, procs, mat, n);
    // gather, send and recieve are the same buffer for rank 0 so skip it with MPI_IN_PLACE
    MPI_Gather(mat.data(), (n*n) / procs, MPI_FLOAT, mat.data(), (n*n) / procs, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // print result
    if (rank == 0) {
        timing = MPI_Wtime() - timing;
        cout << "Ending timer..." << endl;
        MatrixXf L;
        MatrixXf U;
        split_lu(mat, L, U);

        // check correctness
        MatrixXf reread = ext == "csv" ? load_csv(filename) : read_matrix_market(filename);
        bool correct = reread.isApprox((L*U));
        cout << "Correct answer: " << (correct ? "Yes" : "No") << endl;
        cout << "Time: " << timing << endl;
    }
    MPI_Finalize();
    return 0;
}