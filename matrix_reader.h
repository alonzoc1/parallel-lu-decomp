#ifndef MATRIX_READER_H
#define MATRIX_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "fmm/app/Eigen.hpp"

int** readCSV(std::string& filename, int size);

Eigen::MatrixXf read_matrix_market(std::string filename);

#endif