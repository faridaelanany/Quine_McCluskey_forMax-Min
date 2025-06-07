#include "qm.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <map>
#include <set>
#include <vector>
#include <bitset>

using namespace std;

// Initializes the QM minimizer with number of variables (up to 20)
QM::QM(int variables) : VARIABLES(variables) {
    if (variables < 1 || variables > 20) {
        throw invalid_argument("Number of variables must be between 1 and 20.");
    }
}

// Converts a decimal number to binary string representation
string QM::decToBin(int n) {
    if (VARIABLES == 0) return "";
    // Uses bitset to handle conversion and pads to 20 bits, then truncates
    return bitset<20>(n).to_string().substr(20 - VARIABLES);
}

//Checks if two terms differ by exactly one bit
bool QM::isGreyCode(const string& a, const string& b) {
    if (a.length() != b.length()) return false;
    int diffCount = 0;
    for (size_t i = 0; i < a.length(); i++) {
        if (a[i] != b[i]) {
            diffCount++;
            if (diffCount > 1) return false; // Early exit if more than 1 difference
        }
    }
    return diffCount == 1;
}

// Combines two terms by replacing differing bit with '-'
string QM::combineTerms(const string& a, const string& b) {
    string result;
    for (size_t i = 0; i < a.length(); i++) {
        result += (a[i] == b[i]) ? a[i] : '-'; // Keep matching bits, else insert don't-care
    }
    return result;
}

// Checks if a term covers a specific binary minterm
bool QM::covers(const string& term, const string& binaryMinterm) {
    if (term.length() != binaryMinterm.length()) return false;
    for (size_t i = 0; i < term.length(); i++) {
        // For each bit: if not don't-care and doesn't match, return false
        if (term[i] != '-' && term[i] != binaryMinterm[i]) {
            return false;
        }
    }
    return true;
}

// Checks if a term covers a decimal minterm
bool QM::covers(const string& term, int minterm) {
    return covers(term, decToBin(minterm));
}

// Converts maxterms to minterms using complementation
vector<int> QM::convertMaxtermsToMinterms(const vector<int>& maxterms) {
    vector<int> minterms;
    int totalTerms = (1 << VARIABLES);
    set<int> maxtermSet(maxterms.begin(), maxterms.end());

    // All terms not in maxterm list are minterms
    for (int i = 0; i < totalTerms; i++) {
        if (!maxtermSet.count(i)) {
            minterms.push_back(i);
        }
    }
    return minterms;
}

// our function to generate all prime implicants
void QM::generatePrimeImplicants() {
    // Combine minterms and don't-cares, remove duplicates
    vector<int> allTerms = mintermList;
    allTerms.insert(allTerms.end(), dontCareList.begin(), dontCareList.end());
    sort(allTerms.begin(), allTerms.end());
    allTerms.erase(unique(allTerms.begin(), allTerms.end()), allTerms.end());

    // Convert all terms to binary strings
    vector<string> binaryTerms;
    for (int term : allTerms) {
        binaryTerms.push_back(decToBin(term));
    }

    if (binaryTerms.empty()) {
        primeImplicants.clear();
        return;
    }

    // Group terms by number of 1s (key = count of 1s, value = list of terms)
    map<int, vector<string>> groups;
    for (const string& term : binaryTerms) {
        int oneCount = count(term.begin(), term.end(), '1');
        groups[oneCount].push_back(term);
    }

    // Track which terms have been combined
    map<string, bool> isCombined;
    for (const string& term : binaryTerms) {
        isCombined[term] = false;
    }

    vector<string> currentPIs = binaryTerms;
    bool changed = true;

    // Main combining loop
    while (changed) {
        changed = false;
        vector<string> nextPIs;
        set<string> marked; // Terms that get combined

        // Compare adjacent groups (terms differing by one 1 count)
        for (auto it = groups.begin(); it != groups.end(); ++it) {
            auto next = it;
            ++next;
            if (next == groups.end()) break;

            // Compare all terms in current group with next group
            for (const string& term1 : it->second) {
                for (const string& term2 : next->second) {
                    if (isGreyCode(term1, term2)) {
                        string combined = combineTerms(term1, term2);
                        if (find(nextPIs.begin(), nextPIs.end(), combined) == nextPIs.end()) {
                            nextPIs.push_back(combined);
                        }
                        marked.insert(term1);
                        marked.insert(term2);
                        changed = true;
                    }
                }
            }
        }

        // Add unmarked terms to prime implicants (they couldn't be combined further)
        for (const string& term : currentPIs) {
            if (marked.find(term) == marked.end()) {
                if (find(primeImplicants.begin(), primeImplicants.end(), term) == primeImplicants.end()) {
                    primeImplicants.push_back(term);
                }
            }
        }

        // Prepare for next iteration
        if (changed) {
            groups.clear();
            currentPIs = nextPIs;
            // Regroup the new combined terms
            for (const string& term : nextPIs) {
                int oneCount = count(term.begin(), term.end(), '1');
                groups[oneCount].push_back(term);
            }
        }
    }

    // If no combining occurred, each term is its own PI
    if (primeImplicants.empty()) {
        primeImplicants = binaryTerms;
    }

    // Remove duplicate prime implicants
    sort(primeImplicants.begin(), primeImplicants.end());
    primeImplicants.erase(unique(primeImplicants.begin(), primeImplicants.end()), primeImplicants.end());

    // Build coverage map (only for minterms, not don't-cares)
    implicantCoverage.clear();
    for (const string& pi : primeImplicants) {
        set<int> covered;
        for (int minterm : mintermList) {
            if (covers(pi, minterm)) {
                covered.insert(minterm);
            }
        }
        if (!covered.empty()) {
            implicantCoverage[pi] = covered;
        }
    }
}

//Identifies essential prime implicants
void QM::findEssentialPrimeImplicants() {
    if (primeImplicants.empty()) {
        essentialPrimeImplicants.clear();
        minimalSolutions.clear();
        return;
    }

    // Create coverage map: minterm → list of PIs that cover it
    map<int, vector<string>> mintermCoverage;
    for (const auto& pi : primeImplicants) {
        for (int minterm : mintermList) {
            if (covers(pi, minterm)) {
                mintermCoverage[minterm].push_back(pi);
            }
        }
    }

    // Find essential PIs (terms that are the only cover for some minterm)
    essentialPrimeImplicants.clear();
    set<string> essentialPIs;
    set<int> coveredMinterms;

    for (const auto& entry : mintermCoverage) {
        if (entry.second.size() == 1) { // Only one PI covers this minterm → essential
            const string& essentialPI = entry.second[0];
            if (essentialPIs.insert(essentialPI).second) {
                essentialPrimeImplicants.push_back(essentialPI);
                // Mark all minterms this essential PI covers
                for (int m : implicantCoverage[essentialPI]) {
                    coveredMinterms.insert(m);
                }
            }
        }
    }

    // Find minterms not covered by essential PIs
    set<int> uncoveredMinterms;
    for (int m : mintermList) {
        if (coveredMinterms.find(m) == coveredMinterms.end()) {
            uncoveredMinterms.insert(m);
        }
    }

    uncoveredMintermsAfterEPI = vector<int>(uncoveredMinterms.begin(), uncoveredMinterms.end());

    if (uncoveredMinterms.empty()) {
        minimalSolutions.clear();
        return;
    }

    // Find remaining PIs (non-essential ones that cover uncovered minterms)
    map<string, set<int>> remainingCoverage;
    vector<string> remainingPIs;
    for (const string& pi : primeImplicants) {
        if (essentialPIs.find(pi) != essentialPIs.end()) continue;

        set<int> coverage;
        for (int m : uncoveredMinterms) {
            if (covers(pi, m)) {
                coverage.insert(m);
            }
        }
        if (!coverage.empty()) {
            remainingPIs.push_back(pi);
            remainingCoverage[pi] = coverage;
        }
    }

    if (remainingPIs.empty()) {
        minimalSolutions.clear();
        return;
    }

    // Use Petrick's method to select minimal set of remaining PIs
    petricksMethod(remainingPIs, remainingCoverage, uncoveredMinterms);
}

/* Petrick's method for selecting minimal cover of remaining minterms */
void QM::petricksMethod(const vector<string>& remainingPIs,
                       const map<string, set<int>>& remainingCoverage,
                       const set<int>& uncoveredMinterms) {
    minimalSolutions.clear();

    if (remainingPIs.empty() || uncoveredMinterms.empty()) {
        return;
    }

    // Create product-of-sums: for each minterm, list of PIs that cover it
    map<int, vector<string>> mintermToPIs;
    for (int m : uncoveredMinterms) {
        for (const string& pi : remainingPIs) {
            if (remainingCoverage.at(pi).count(m)) {
                mintermToPIs[m].push_back(pi);
            }
        }
    }

    // Initialize solutions with first minterm's PIs
    vector<vector<string>> solutions;
    if (!mintermToPIs.empty()) {
        for (const string& pi : mintermToPIs.begin()->second) {
            solutions.push_back({pi});
        }
    }

    // Multiply solutions (AND operation between product terms)
    for (auto it = next(mintermToPIs.begin()); it != mintermToPIs.end(); ++it) {
        vector<vector<string>> newSolutions;
        for (const vector<string>& sol : solutions) {
            for (const string& pi : it->second) {
                vector<string> newSol = sol;
                if (find(newSol.begin(), newSol.end(), pi) == newSol.end()) {
                    newSol.push_back(pi);
                }
                // Remove duplicates and keep sorted
                sort(newSol.begin(), newSol.end());
                newSol.erase(unique(newSol.begin(), newSol.end()), newSol.end());
                // Add if not already present
                if (find(newSolutions.begin(), newSolutions.end(), newSol) == newSolutions.end()) {
                    newSolutions.push_back(newSol);
                }
            }
        }
        solutions = newSolutions;

        // Early exit if we find a solution with just one PI
        bool hasMinimal = false;
        for (const auto& sol : solutions) {
            if (sol.size() == 1) {
                hasMinimal = true;
                break;
            }
        }
        if (hasMinimal) break;
    }

    // Find all solutions with minimal size
    if (!solutions.empty()) {
        size_t minSize = solutions[0].size();
        for (const auto& sol : solutions) {
            if (sol.size() < minSize) {
                minSize = sol.size();
            }
        }

        for (const auto& sol : solutions) {
            if (sol.size() == minSize) {
                minimalSolutions.push_back(sol);
            }
        }
    }
}

// Converts binary representation to Boolean expression
string QM::binaryToExpression(const string& binary) {
    if (VARIABLES == 0) return "";

    vector<string> vars;
    for (int i = 0; i < VARIABLES; i++) {
        vars.push_back(string(1, 'A' + i)); // Variables are A, B, C, etc.
    }

    string expression;
    for (size_t i = 0; i < binary.length(); i++) {
        if (binary[i] == '0') {
            expression += vars[i] + "'"; // Complemented variable
        }
        else if (binary[i] == '1') {
            expression += vars[i];       // Normal variable
        }
        // '-' is ignored (don't-care)
    }

    return expression;
}

// Prints a table showing coverage of each prime implicant
void QM::printCoverageTable() {
    cout << "\nPrime Implicants Coverage Table:\n";
    cout << "| Prime Implicant | Binary Representation | Covers Minterms | Covers Don't-cares |\n";
    cout << "|-----------------|-----------------------|-----------------|--------------------|\n";

    for (const string& pi : primeImplicants) {
        // Find covered minterms
        set<int> coveredMinterms;
        for (int m : mintermList) {
            if (covers(pi, m)) coveredMinterms.insert(m);
        }

        // Find covered don't-cares
        set<int> coveredDontCares;
        for (int dc : dontCareList) {
            if (covers(pi, dc)) coveredDontCares.insert(dc);
        }

        // Format table row
        cout << "| " << setw(15) << binaryToExpression(pi) << " | "
            << setw(21) << pi << " | ";

        if (!coveredMinterms.empty()) {
            cout << "{";
            bool first = true;
            for (int m : coveredMinterms) {
                if (!first) cout << ", ";
                cout << m;
                first = false;
            }
            cout << "}";
        }
        else {
            cout << "None";
        }

        cout << " | ";

        if (!coveredDontCares.empty()) {
            cout << "{";
            bool first = true;
            for (int dc : coveredDontCares) {
                if (!first) cout << ", ";
                cout << dc;
                first = false;
            }
            cout << "}";
        }
        else {
            cout << "None";
        }

        cout << " |\n";
    }
}

// Parses from input string ( they are comma seperated)
vector<int> QM::parseIntegers(const string& input) {
    vector<int> result;
    istringstream iss(input);
    string token;
    while (getline(iss, token, ',')) {
        try {
            // Remove whitespace
            token.erase(remove_if(token.begin(), token.end(), ::isspace), token.end());
            if (!token.empty()) {
                int value = stoi(token);
                result.push_back(value);
            }
        }
        catch (const invalid_argument& e) {
            cerr << "Warning: Invalid term '" << token << "' will be ignored.\n";
        }
        catch (const out_of_range& e) {
            cerr << "Warning: Term '" << token << "' is out of range and will be ignored.\n";
        }
    }
    return result;
}

// Validates input minterms and don't-cares
bool QM::validateInput() {
    // Remove duplicates
    sort(mintermList.begin(), mintermList.end());
    mintermList.erase(unique(mintermList.begin(), mintermList.end()), mintermList.end());

    sort(dontCareList.begin(), dontCareList.end());
    dontCareList.erase(unique(dontCareList.begin(), dontCareList.end()), dontCareList.end());

    // Validate term ranges
    int maxTerm = (1 << VARIABLES) - 1;
    bool valid = true;

    // Check for overlapping minterms and don't-cares
    vector<int> intersection;
    set_intersection(mintermList.begin(), mintermList.end(),
                    dontCareList.begin(), dontCareList.end(),
                    back_inserter(intersection));

    if (!intersection.empty()) {
        cerr << "Error: The following terms appear in both minterm and don't-care lists: ";
        for (size_t i = 0; i < intersection.size(); i++) {
            if (i != 0) cerr << ", ";
            cerr << intersection[i];
        }
        cerr << endl;
        valid = false;
    }

    // Check minterms are in valid range
    vector<int> validMterms;
    for (int m : mintermList) {
        if (m < 0 || m > maxTerm) {
            cerr << "Error: Minterm " << m << " is out of range (0-" << maxTerm << ")\n";
            valid = false;
        } else {
            validMterms.push_back(m);
        }
    }
    mintermList = validMterms;

    // Check don't-cares are in valid range
    vector<int> validDCs;
    for (int dc : dontCareList) {
        if (dc < 0 || dc > maxTerm) {
            cerr << "Error: Don't-care term " << dc << " is out of range (0-" << maxTerm << ")\n";
            valid = false;
        } else {
            validDCs.push_back(dc);
        }
    }
    dontCareList = validDCs;

    return valid;
}

// Generates Verilog module implementing the minimized function
void QM::printVerilogModule() {
    cout << "\nVerilog Module (Structural):\n";
    // Module declaration
    cout << "module minimized_function(";
    for (int i = 0; i < VARIABLES; i++) {
        if (i != 0) cout << ", ";
        cout << char('A' + i);
    }
    cout << ", F);\n";

    // Input/output declarations
    cout << "  input ";
    for (int i = 0; i < VARIABLES; i++) {
        if (i != 0) cout << ", ";
        cout << char('A' + i);
    }
    cout << ";\n";
    cout << "  output F;\n\n";

    if (essentialPrimeImplicants.empty()) {
        // Handle constant outputs
        if (mintermList.empty()) {
            cout << "  // Constant 0 output\n";
            cout << "  buf(F, 1'b0);\n";
        } else {
            cout << "  // Constant 1 output\n";
            cout << "  buf(F, 1'b1);\n";
        }
    } else {
        // Declare wires for intermediate signals
        for (size_t i = 0; i < essentialPrimeImplicants.size(); i++) {
            cout << "  wire p" << i << ";\n";
        }
        if (!minimalSolutions.empty()) {
            for (size_t i = 0; i < minimalSolutions[0].size(); i++) {
                cout << "  wire s" << i << ";\n";
            }
        }
        cout << "  wire or_out;\n\n";

        // Generate NOT gates for complemented inputs
        for (size_t i = 0; i < essentialPrimeImplicants.size(); i++) {
            const string& pi = essentialPrimeImplicants[i];
            for (size_t j = 0; j < pi.length(); j++) {
                if (pi[j] == '0') {
                    cout << "  not not_" << char('a' + j) << "_p" << i << "(not_" << char('a' + j) << "_p" << i << ", " << char('A' + j) << ");\n";
                }
            }
        }
        if (!minimalSolutions.empty()) {
            for (size_t i = 0; i < minimalSolutions[0].size(); i++) {
                const string& pi = minimalSolutions[0][i];
                for (size_t j = 0; j < pi.length(); j++) {
                    if (pi[j] == '0') {
                        cout << "  not not_" << char('a' + j) << "_s" << i << "(not_" << char('a' + j) << "_s" << i << ", " << char('A' + j) << ");\n";
                    }
                }
            }
        }
        cout << "\n";

        // Generate AND gates for product terms (essential PIs)
        for (size_t i = 0; i < essentialPrimeImplicants.size(); i++) {
            const string& pi = essentialPrimeImplicants[i];
            cout << "  and and_p" << i << "(p" << i;

            for (size_t j = 0; j < pi.length(); j++) {
                if (pi[j] == '0') {
                    cout << ", not_" << char('a' + j) << "_p" << i;
                } else if (pi[j] == '1') {
                    cout << ", " << char('A' + j);
                }
            }
            cout << ");\n";
        }

        // Generate AND gates for secondary PIs (if any)
        if (!minimalSolutions.empty()) {
            for (size_t i = 0; i < minimalSolutions[0].size(); i++) {
                const string& pi = minimalSolutions[0][i];
                cout << "  and and_s" << i << "(s" << i;

                for (size_t j = 0; j < pi.length(); j++) {
                    if (pi[j] == '0') {
                        cout << ", not_" << char('a' + j) << "_s" << i;
                    } else if (pi[j] == '1') {
                        cout << ", " << char('A' + j);
                    }
                }
                cout << ");\n";
            }
        }
        cout << "\n";

        // Generate OR gate combining all product terms
        cout << "  or or_gate(or_out";
        for (size_t i = 0; i < essentialPrimeImplicants.size(); i++) {
            cout << ", p" << i;
        }
        if (!minimalSolutions.empty()) {
            for (size_t i = 0; i < minimalSolutions[0].size(); i++) {
                cout << ", s" << i;
            }
        }
        cout << ");\n";

        // Final output assignment
        cout << "  buf(F, or_out);\n";
    }

    cout << "endmodule\n";
}

// our minimization function that coordinates all steps ( output function)
void QM::minimize() {
    if (!validateInput()) {
        cerr << "Input validation failed. Cannot proceed with minimization.\n";
        return;
    }

    if (mintermList.empty() && dontCareList.empty()) {
        cout << "No minterms or don't-care terms provided. Nothing to minimize.\n";
        return;
    }

    // Execute minimization steps
    generatePrimeImplicants();
    findEssentialPrimeImplicants();

    // Print results
    cout << "\n--- Quine-McCluskey Minimization Results ---\n";
    cout << "Number of variables: " << VARIABLES << "\n";
    cout << "Minterms: ";
    for (size_t i = 0; i < mintermList.size(); i++) {
        if (i != 0) cout << ", ";
        cout << mintermList[i];
    }
    if (mintermList.empty()) cout << "None";
    cout << "\nDon't-care terms: ";
    for (size_t i = 0; i < dontCareList.size(); i++) {
        if (i != 0) cout << ", ";
        cout << dontCareList[i];
    }
    if (dontCareList.empty()) cout << "None";
    cout << "\n";

    // Print all prime implicants
    cout << "\nAll Prime Implicants (" << primeImplicants.size() << "):\n";
    for (const string& pi : primeImplicants) {
        cout << binaryToExpression(pi) << " (" << pi << ")\n";
    }

    // Print essential prime implicants
    cout << "\nEssential Prime Implicants (" << essentialPrimeImplicants.size() << "):\n";
    for (const string& epi : essentialPrimeImplicants) {
        cout << binaryToExpression(epi) << " (" << epi << ")\n";
    }

    // Print uncovered minterms (if any)
    if (!uncoveredMintermsAfterEPI.empty()) {
        cout << "\nMinterms not covered by essential PIs: ";
        for (size_t i = 0; i < uncoveredMintermsAfterEPI.size(); i++) {
            if (i != 0) cout << ", ";
            cout << uncoveredMintermsAfterEPI[i];
        }
        cout << "\n";
    }

    // Print coverage table
    printCoverageTable();

    // Print minimized expression
    cout << "\nMinimized Boolean Expression: ";
    if (essentialPrimeImplicants.empty()) {
        if (mintermList.empty()) {
            cout << "0 (No minterms)";
        } else {
            cout << "1 (All minterms covered by don't-cares)";
        }
    } else {
        bool first = true;
        // Print essential PIs
        for (const string& epi : essentialPrimeImplicants) {
            if (!first) cout << " + ";
            cout << binaryToExpression(epi);
            first = false;
        }

        // Print additional PIs from minimal solution (if any)
        if (!minimalSolutions.empty()) {
            for (const string& pi : minimalSolutions[0]) {
                if (!first) cout << " + ";
                cout << binaryToExpression(pi);
                first = false;
            }
        }
    }

    // Print alternative solutions (if any)
    if (!minimalSolutions.empty() && minimalSolutions.size() > 1) {
        cout << "\n\nAlternative minimal solutions (" << minimalSolutions.size() << "):\n";
        for (size_t i = 0; i < minimalSolutions.size(); i++) {
            cout << "Solution " << (i+1) << ": ";
            bool firstTerm = true;
            for (const string& pi : essentialPrimeImplicants) {
                if (!firstTerm) cout << " + ";
                cout << binaryToExpression(pi);
                firstTerm = false;
            }
            for (const string& pi : minimalSolutions[i]) {
                if (!firstTerm) cout << " + ";
                cout << binaryToExpression(pi);
                firstTerm = false;
            }
            cout << "\n";
        }
    }

    cout << endl;

    // Generate Verilog implementation
    printVerilogModule();
}

// Reads minimization problem from file
void QM::readFromFile(const string& filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    string line;
    int lineNum = 0;
    bool isMaxtermFile = false;

    // Reset previous data
    mintermList.clear();
    dontCareList.clear();
    primeImplicants.clear();
    essentialPrimeImplicants.clear();
    implicantCoverage.clear();
    minimalSolutions.clear();
    uncoveredMintermsAfterEPI.clear();

    // Read number of variables (first line)
    if (getline(infile, line)) {
        lineNum++;
        try {
            VARIABLES = stoi(line);
            if (VARIABLES < 1 || VARIABLES > 20) {
                throw out_of_range("Number of variables must be between 1 and 20");
            }
        }
        catch (const exception& e) {
            throw runtime_error("Line " + to_string(lineNum) + ": Invalid number of variables. " + e.what());
        }
    }
    else {
        throw runtime_error("File is empty");
    }

    // Read second line to determine if minterms or maxterms
    if (getline(infile, line)) {
        lineNum++;
        string lowercaseLine = line;
        transform(lowercaseLine.begin(), lowercaseLine.end(), lowercaseLine.begin(), ::tolower);

        if (lowercaseLine.find("maxterms") != string::npos) {
            isMaxtermFile = true;
            // Read maxterms from next line
            if (getline(infile, line)) {
                lineNum++;
                vector<int> maxterms = parseIntegers(line);
                mintermList = convertMaxtermsToMinterms(maxterms);
            }
            else {
                throw runtime_error("Missing maxterms line");
            }
        }
        else if (lowercaseLine.find("minterms") != string::npos) {
            // Read minterms from next line
            if (getline(infile, line)) {
                lineNum++;
                mintermList = parseIntegers(line);
            }
            else {
                throw runtime_error("Missing minterms line");
            }
        }
        else {
            // Assume line contains minterms without keyword
            mintermList = parseIntegers(line);
        }
    }
    else {
        throw runtime_error("Missing minterms/maxterms line");
    }

    // Read don't-care terms if present (third line)
    if (getline(infile, line)) {
        lineNum++;
        dontCareList = parseIntegers(line);
    }

    infile.close();

    // Remove duplicates
    sort(mintermList.begin(), mintermList.end());
    mintermList.erase(unique(mintermList.begin(), mintermList.end()), mintermList.end());

    sort(dontCareList.begin(), dontCareList.end());
    dontCareList.erase(unique(dontCareList.begin(), dontCareList.end()), dontCareList.end());

    if (isMaxtermFile) {
        cout << "Processed maxterm file. Converted maxterms to minterms for minimization.\n";
        // Remove don't-care terms from minterm list if they overlap
        vector<int> intersection;
        set_intersection(mintermList.begin(), mintermList.end(),
                        dontCareList.begin(), dontCareList.end(),
                        back_inserter(intersection));

        if (!intersection.empty()) {
            vector<int> newMinterms;
            set_difference(mintermList.begin(), mintermList.end(),
                         dontCareList.begin(), dontCareList.end(),
                         back_inserter(newMinterms));
            mintermList = newMinterms;

            cout << "Note: Removed " << intersection.size()
                 << " don't-care terms from minterm list.\n";
        }
    }
}