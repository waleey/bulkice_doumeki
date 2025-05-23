#include "OMSimUtils.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <algorithm>

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
        ElementData element;
        std::string field;

        // Read tab-separated values correctly
        std::getline(ss, element.isotopeName, '\t');  // Isotope name
        std::getline(ss, field, '\t'); element.atomicNumber = std::stoi(field); // Atomic Number
        std::getline(ss, field, '\t'); element.massNumber = std::stoi(field);   // Mass Number
        std::getline(ss, field, '\t'); element.excitationEnergy = std::stod(field); // Excitation Energy
        std::getline(ss, field, '\t'); element.totalAngularMomentum = std::stoi(field); // Total Angular Momentum
        std::getline(ss, field, '\t'); element.branchingRatio = std::stod(field); // Branching Ratio

        if (ss.fail()) {
            std::cerr << "Error: Failed to parse line: " << line << std::endl;
        } else {
            data.push_back(element);
        }
    }

    file.close();
    return data;
}

std::string sanitize_for_filename(double value) {
    std::string formatted = boost::str(boost::format("%.0E") % value); // Scientific, uppercase E, 1 decimal
    // Remove the decimal point (e.g., 5.0E-01 → 50E-01)
    formatted.erase(std::remove(formatted.begin(), formatted.end(), '.'), formatted.end());
    return formatted;
}