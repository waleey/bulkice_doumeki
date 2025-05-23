#ifndef OMSimUtils_h
#define OMSimUtils_h

#include <string>
#include <vector>

// Define the ElementData struct
struct ElementData {
    std::string isotopeName;  // E.g., "U238"
    int atomicNumber;  // E.g., 92
    int massNumber;    // E.g., 238
    double excitationEnergy;   // in [keV] E.g., 0.0, 73.92, ...
    int totalAngularMomentum;     // in [1/2] E.g., 0, 3, 6
    double branchingRatio; // E.g., 1.0
};

// Declare the loadDataFromFile function
std::vector<ElementData> loadDataFromFile(const std::string& filename);
std::string sanitize_for_filename(double value);

#endif // OmSimUtils_h