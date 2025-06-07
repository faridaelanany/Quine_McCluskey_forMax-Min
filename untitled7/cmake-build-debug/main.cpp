#include "qm.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;

int main() {
    try {
        string filename;
        cout << "Quine-McCluskey Boolean Function Minimizer\n";
        cout << "Supports functions with up to 20 variables\n";
        cout << "Enter input file name: ";
        getline(cin, filename);

        ifstream testFile(filename);
        if (!testFile.is_open()) {
            cerr << "Error: Could not open file " << filename << endl;
            return 1;
        }

        string firstLine;
        getline(testFile, firstLine);
        testFile.close();

        int numVariables;
        try {
            numVariables = stoi(firstLine);
            if (numVariables < 1 || numVariables > 20) {
                throw out_of_range("Number of variables must be between 1 and 20");
            }
        }
        catch (const exception& e) {
            cerr << "Error in input file: " << e.what() << endl;
            cerr << "First line must be number of variables (1-20)" << endl;
            return 1;
        }

        QM qm(numVariables);
        try {
            qm.readFromFile(filename);
        }
        catch (const exception& e) {
            cerr << "Error reading file: " << e.what() << endl;
            cerr << "File format must be:\n";
            cerr << "Line 1: Number of variables\n";
            cerr << "Line 2: 'maxterms' (optional) followed by terms\n";
            cerr << "Line 3: Don't-care terms (optional)\n";
            return 1;
        }

        qm.minimize();
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}