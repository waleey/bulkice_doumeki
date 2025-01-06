#include "OMSimUtils.hh"
#include <iostream>
#include <fstream>
#include <sstream>

// Implement the loadDataFromFile function
std::vector<ElementData> loadDataFromFile(const std::string& filename) {
    std::vector<ElementData> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    std::string line;
    bool isHeader = true; // Assume the first line is a header

    while (std::getline(file, line)) {
        // Skip header
        if (isHeader) {
            isHeader = false;
            continue;
        }

        std::stringstream ss(line);
        std::string token;
        ElementData element;

        // Parse the line assuming comma or tab separation
        std::getline(ss, element.name, '\t');               // First column: Name

        if (ss.fail()) {
            std::cerr << "Failed to parse with '\t' delimiter. ss state: " << ss.str() << std::endl;
        }
        
        ss >> element.atomicNumber; ss.ignore(1, ',');     // Second column: Atomic Number
        ss >> element.massNumber; ss.ignore(1, ',');       // Third column: Mass Number
        ss >> element.excitationEnergy; ss.ignore(1, ',');         // Fourth column: Half-Life
        ss >> element.totalAngularMomentum; ss.ignore(1, ',');        // Fifth column: Decay Mode
        ss >> element.branchingRatio;                      // Sixth column: Branching Ratio

        if (!ss.fail()) {
            data.push_back(element);
        } else {
            std::cerr << "Error: Failed to parse line: " << line << std::endl;
        }
    }

    file.close();
    return data;
}
