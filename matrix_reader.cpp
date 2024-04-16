
#include "matrix_reader.h"

using namespace std;

int** vectorTo2DArray(std::vector<std::vector<int>>& data, int row_count, int col_count) {
    int** result = new int*[row_count];
    for (int i = 0; i < row_count; i++) {
        result[i] = new int[col_count];
        for (int j = 0; j < col_count; j++) {
            result[i][j] = data[i][j];
        }
    }
    return result;
}

std::vector<int> split_row(const string& string_stream, char delim) {
    std::vector<int> toks;
    string tok;
    istringstream tokenStream(string_stream);
    while (getline(tokenStream, tok, delim)) {
        toks.push_back(stoi(tok));
    }
    return toks;
}

int** readCSV(string& filename, int size) {
    ifstream file(filename);
    // We use vectors to dynamically read in the data
    std::vector<std::vector<int>> data;
    if (file) {
        string line;
        while (getline(file, line)) {
            std::vector<int> row = split_row(line, ',');
            data.push_back(row);
        }
        file.close();
    } else {
        cerr << "Unable to read file " << filename << endl;
    }
    return vectorTo2DArray(data, size, size);
}

/* // Uncomment this if you want to test matrix_reader on its own
int main (int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "To use this function, provide a <filename.csv> and <matrix_size>" << endl;
        cerr << "Usage: " << argv[0] << "<filename.csv>" << " " << "<matrix_size>" << endl;
        return 1;
    }
    string filename = argv[1];
    int** result = readCSV(filename, stoi(argv[2]));
    for (int i = 0; i < stoi(argv[2]); i++) {
        delete[] result[i];
    }
    delete[] result;
    return 0;
}
*/