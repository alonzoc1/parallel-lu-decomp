# parallel-lu-decomp
Our parallel implementation of Guassian LU decomposition for CS605

# Build instructions
To build gaussian serial implementation, use `g++ gaussian_serial.cpp matrix_reader.cpp -o gaussian_serial`. Run it with `./gaussian_serial filename.csv` or `./gaussian_serial filename.mtx`.

To build MPI parallel gaussian implementation, use `mpic++ gaussian_parallel_mpi.cpp matrix_reader.cpp -o gaussian_mpi`. Run it with `mpiexec -np NUM_PROCS ./gaussian_mpi filename.csv MAT_SIZE` or `mpiexec -np NUM_PROCS ./gaussian_mpi filename.mtx MAT_SIZE` where `NUM_PROCS` is the number of processes to provide the algorithm, and MAT_SIZE is the size of the matrix, which needs to be provided beforehand. This is because this parallel implementation has to allocate all the contiguous memory before loading the matrix.

Eigen and Fast Matrix Market headers are included in this repository. For serial timing, boost library is required (commonly found, search for steps specific to your OS, or see generic steps here: https://www.boost.org/doc/libs/1_55_0/doc/html/bbv2/installation.html)

Fast Matrix Market may require certain dependencies to be installed on your system (https://github.com/alugowski/fast_matrix_market)
