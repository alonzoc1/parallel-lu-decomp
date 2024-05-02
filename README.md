# parallel-lu-decomp
Our parallel implementation of Guassian LU decomposition for CS605

# Build instructions
To build gaussian serial implementation, use `g++ gaussian_serial.cpp matrix_reader.cpp -o gaussian_serial`. Run it with `./gaussian_serial filename.csv` or `./gaussian_serial filename.mtx`.

To build MPI parallel gaussian implementation, use `mpic++ gaussian_parallel_mpi.cpp matrix_reader.cpp -o gaussian_mpi`. Run it with `mpiexec -np NUM_PROCS ./gaussian_mpi filename.csv MAT_SIZE` or `mpiexec -np NUM_PROCS ./gaussian_mpi filename.mtx MAT_SIZE` where `NUM_PROCS` is the number of processes to provide the algorithm, and MAT_SIZE is the size of the matrix, which needs to be provided beforehand. This is because this parallel implementation has to allocate all the contiguous memory before loading the matrix.

Finally, a required component to run our automated testing code is the matrix generator. Please compile it before hand using `g++ matrix_gen.cpp -o gen`. It is used to generate random matrices.

Eigen and Fast Matrix Market headers are included in this repository. For serial timing, boost library is required (commonly found, search for steps specific to your OS, or see generic steps here: https://www.boost.org/doc/libs/1_55_0/doc/html/bbv2/installation.html)

Fast Matrix Market may require certain dependencies to be installed on your system (https://github.com/alugowski/fast_matrix_market)

# Test cases and how to replicate tests
## MPI and serial
Results are stored in raw_serial.txt and raw_mpi.txt at this point. `run_tests.cpp` can be compiled (only uses boost timer) and run with `./run_tests TEST_TYPE AVG_OVER` where TEST_TYPE is `serial` or `mpi`, and AVG_OVER is how many iterations per tested size. Tested matrix sizes are 64, 256, 512, 1024, 2048. For each iteration a random matrix of that size is generated using `matrix_gen.cpp` (expected to be compiled as `gen`), and is run using the compiled serial or mpi implementation.

For example, to test each size 3 times for serial, run `./run_tests serial 3`, which will output to `raw_serial.txt` (appends to file, so delete or rename this file for a clean run).

When MPI is run in this way, automatically tests with process counts 2, 4, and 8.
