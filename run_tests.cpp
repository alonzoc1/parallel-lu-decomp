
#include <cstdlib>
#include <iostream>
#ifndef BOOST_TIMER_ENABLE_DEPRECATED
#define BOOST_TIMER_ENABLE_DEPRECATED
#endif
#include <boost/timer.hpp>

using namespace std;

int AVG_OVER;

string build_test_data_filename(int iteration, int size) {
    string test_data_prefix = "./data/test_data_";
    test_data_prefix += to_string(size) + "_" + to_string(iteration) + ".csv";
    return test_data_prefix;
}

/* Generates random matrices for sizes 64, 256, 512, 1024, and 2048, then tests them the given amount of times
Saves results to results_TYPE.txt, and raw outputs to raw_TYPE.txt, where TYPE is the test type
*/
int main(int argc, char* argv[]) {
    string usage = "Usage: ./run_tests TEST_TYPE AVG_OVER\nWhere TEST_TYPE is serial or mpi, and AVG_OVER is how many iterations to run each matrix";
    if (argc != 3) {
        cout << usage << endl;
        return 1;
    }
    string test_type = argv[1];
    AVG_OVER = stoi(argv[2]);
    string raw_filename = "raw_" + test_type + ".txt";
    string results_filename = "results_" + test_type + ".txt";
    int test_sizes[4] = {64, 256, 512, 1024};
    // Generate data first
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < AVG_OVER; j++) {
            cout << "Generating data for " << test_sizes[i] << " iteration " << j << endl;
            if (system(("./gen " + build_test_data_filename(j, test_sizes[i]) + " " + to_string(test_sizes[i])).c_str()) != 0) {
                cout << "Generating data failed, exiting..." << endl;
                return 1;
            }
            //sleep(1); // let seed increment
        }
    }
    // Run tests
    for (int i = 0; i < 4; i++) {
        int test_size = test_sizes[i];
        // Start timer
        cout << "Starting timer for size: " << test_size << endl;
        boost::timer myTimer;
        for (int iter_tracker = 0; iter_tracker < AVG_OVER; iter_tracker++) {
            cout << "--Iteration " << iter_tracker << endl;
            if (test_type == "serial") {
                system(("echo SIZE" + to_string(test_size) + "ITER" + to_string(iter_tracker) + " >> " + raw_filename).c_str());
                if (system(("./gaussian_serial " + build_test_data_filename(iter_tracker, test_size) + " >> " + raw_filename).c_str()) != 0) {
                    cout << "One of the tests failed, exiting..." << endl;
                }
            } else if (test_type == "mpi") {
                int proc_counts[3] = {2, 4, 8};
                for (int proc_i = 0; proc_i < 3; proc_i++) {
                    system(("echo SIZE" + to_string(test_size) + "ITER" + to_string(iter_tracker) + "PROCS" + to_string(proc_counts[proc_i]) + " >> " + raw_filename).c_str());
                    if (system(("mpiexec -np " + to_string(proc_counts[proc_i]) + " -oversubscribe ./gaussian_mpi " + build_test_data_filename(iter_tracker, test_size) + " " + to_string(test_size) + " >> " + raw_filename).c_str()) != 0) {
                        cout << "One of the tests failed, exiting..." << endl;
                    }
                }
            } else if(test_type == "omp"){
                system(("echo SIZE" + to_string(test_size) + "ITER" + to_string(iter_tracker) + " >> " + raw_filename).c_str());
                if (system(("./PCludcmp " + build_test_data_filename(iter_tracker, test_size) + " >> " + raw_filename).c_str()) != 0) {
                    cout << "One of the tests failed, exiting..." << endl;
                }
            }
        }
    }
}