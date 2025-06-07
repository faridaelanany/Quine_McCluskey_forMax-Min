#ifndef QM_H
#define QM_H

#include <vector>
#include <string>
#include <map>
#include <set>

class QM {
public:
    QM(int variables);

    // Main minimization function
    void minimize();

    // File I/O
    void readFromFile(const std::string& filename);

    // Input handling
    std::vector<int> parseIntegers(const std::string& input);
    bool validateInput();

    // Core algorithm functions
    void generatePrimeImplicants();
    void findEssentialPrimeImplicants();
    void applyDominance(std::vector<std::string>& remainingPIs,
                       std::map<std::string, std::set<int>>& remainingCoverage,
                       const std::set<int>& uncoveredMinterms);
    void petricksMethod(const std::vector<std::string>& remainingPIs,
                      const std::map<std::string, std::set<int>>& remainingCoverage,
                      const std::set<int>& uncoveredMinterms);

    // Helper functions
    std::string decToBin(int n);
    bool isGreyCode(const std::string& a, const std::string& b);
    std::string combineTerms(const std::string& a, const std::string& b);
    bool covers(const std::string& term, const std::string& binaryMinterm);
    bool covers(const std::string& term, int minterm);
    std::vector<int> convertMaxtermsToMinterms(const std::vector<int>& maxterms);
    std::string binaryToExpression(const std::string& binary);

    // Output functions
    void printCoverageTable();
    void printVerilogModule();

    // Public member variables for input/output
    std::vector<int> mintermList;
    std::vector<int> dontCareList;
    int VARIABLES;

private:
    std::vector<std::string> primeImplicants;
    std::vector<std::string> essentialPrimeImplicants;
    std::map<std::string, std::set<int>> implicantCoverage;
    std::vector<std::vector<std::string>> minimalSolutions;
    std::vector<int> uncoveredMintermsAfterEPI;
};

#endif // QM_H