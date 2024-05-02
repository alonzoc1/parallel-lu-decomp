
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

float HIGH = 200000.0;
float LOW = -200000.0;
IOFormat CSV(StreamPrecision, DontAlignCols, ", ", "\n");

/* Generate matrix of random values that is always invertible */
void generate_random_matrix(int size, string filename) {
    // set seed from time
    srand((unsigned int) time(0));
    // generate random invertible matrix
    MatrixXf generated = MatrixXf::Random(size, size);
    generated = (generated + MatrixXf::Constant(size,size,1.))*(HIGH-LOW)/2;
    generated = (generated + MatrixXf::Constant(size, size, LOW));
    MatrixXf generated_copy = generated;
    generated_copy.cwiseAbs();
    generated.diagonal() = generated_copy.colwise().sum();
    // save to csv
    ofstream file(filename.c_str());
    file << generated.format(CSV);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        return 1;
    }
    int size = stoi(argv[2]);
    generate_random_matrix(size, argv[1]);
    return 0;
}
