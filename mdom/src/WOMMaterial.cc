#include "WOMMaterial.hh"

/**
*Very rough draft of possible WOM Simulation
*Waly M Z Karim
*8/30/2023
**/

#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

WOMMaterial::WOMMaterial()
{
    GenerateMaterials();
}
WOMMaterial::~WOMMaterial()
{

}
void WOMMaterial::GenerateMaterials()
{
    GetSharedData();
    GenerateGlassMaterial();
    GenerateFillerMaterial();
    GeneratePaintMaterial();
    GenerateTubeMaterial();
    GenerateTubeInsideMaterial();
}
void WOMMaterial::GetSharedData()
{
    quartzDensity = 2.65 * g / cm3; //found from google
    fillerDensity = 1.0 * g / cm3; //for now, used the density of water
    paintThickness = 25 * micrometer;
}
void WOMMaterial::ReadData(const G4String& fileName, std::vector<G4double>& vals)
{
    wLen.clear();
    vals.clear();

    std::fstream inputFile(fileName);

    G4double dumWlen, dumVal;

    if(!inputFile.good())
    {
        std::cerr << "Error reading input file: " << fileName << std::endl
        << "function: ReadData()" << std::endl;
        abort();
    }
    else
    {
        while(inputFile >> dumWlen >> dumVal)
        {
            wLen.push_back(dumWlen);
            vals.push_back(dumVal);
        }
    }
    std::cout << wLen.size() << std::endl;

    inputFile.close();
    inputFile.clear();
}
void WOMMaterial::GenerateGlassMaterial()
{
    //generating quartz

    G4Element* Si = new G4Element("Silicon", "Si", z = 14, a = 28 * g / mole);
    G4Element* O = new G4Element("Oxygen", "O", z = 8, a = 16 * g / mole);

    quartz = new G4Material("Quartz", quartzDensity, nElements = 2);
    quartz -> AddElement(Si, 1);
    quartz -> AddElement(O, 2);

    std::vector<G4double> rIndexEnergy = { 1.775 * eV, 6.212 * eV};
    std::vector<G4double> rIndex = {1.46, 1.46};
    G4double numEntries = rIndexEnergy.size();

    mptQuartz = new G4MaterialPropertiesTable();
    mptQuartz -> AddProperty("RINDEX", rIndexEnergy, rIndex, numEntries);

    quartz -> SetMaterialPropertiesTable(mptQuartz);
}
void WOMMaterial::GenerateFillerMaterial()
{
    //for now, assuming the filler material is Water.
    G4Element* H = new G4Element("Hydrogen", "H", z = 1, a = 1.01 * g / mole);
    G4Element* O = new G4Element("Oxygen", "O", z = 8, a = 16 * g / mole);

    filler = new G4Material("Filler", fillerDensity, nElements = 2);
    filler -> AddElement(H, 2);
    filler -> AddElement(O, 1);

    std::vector<G4double> rIndexEnergy = { 1.775 * eV, 6.212 * eV};
    std::vector<G4double> rIndex = {1.33, 1.33};
    G4double numEntries = rIndexEnergy.size();

    mptFiller = new G4MaterialPropertiesTable();
    mptFiller -> AddProperty("RINDEX", rIndexEnergy, rIndex, numEntries);

    filler -> SetMaterialPropertiesTable(mptFiller);
}
void WOMMaterial::GeneratePaintMaterial()
{
    //generating paint material
    G4Element* Si = new G4Element("Silicon", "Si", z = 14, a = 28 * g / mole);
    G4Element* O = new G4Element("Oxygen", "O", z = 8, a = 16 * g / mole);

    quartzPaint = new G4Material("Quartz", quartzDensity, nElements = 2);
    quartzPaint -> AddElement(Si, 1);
    quartzPaint -> AddElement(O, 2);

    std::vector<G4double> rIndexEnergy = { 1.775 * eV, 6.212 * eV};
    std::vector<G4double> rIndex = {1.46, 1.46};

    std::vector<G4double> absIntensity;
    std::vector<G4double> absLength;
    std::vector<G4double> absEnergy;
    ReadData("/home/waly/bulkice_doumeki/mdom/InputFile/Absorption_WLS.txt", absIntensity);

    for(G4double& val : absIntensity)
    {
        absLength.push_back(- paintThickness / log(val / 100));
    }
    absEnergy = WavelengthToEnergy(wLen);
    absIntensity.clear();

    std::vector<G4double> emIntensity;
    std::vector<G4double> emEnergy;

    ReadData("/home/waly/bulkice_doumeki/mdom/InputFile/Emission_WLS.txt", emIntensity);
    emEnergy = WavelengthToEnergy(wLen);

    mptQuartzPaint = new G4MaterialPropertiesTable();

    mptQuartzPaint -> AddProperty("RINDEX", rIndexEnergy, rIndex, rIndex.size());
    mptQuartzPaint -> AddProperty("WLSABSLENGTH", absEnergy, absLength, absEnergy.size());
    mptQuartzPaint -> AddProperty("WLSCOMPONENT", emEnergy, emIntensity, emIntensity.size());
    mptQuartzPaint -> AddConstProperty("WLSTIMECONSTANT", 1.5 * ns);

    quartzPaint -> SetMaterialPropertiesTable(mptQuartzPaint);
}
void WOMMaterial::GenerateTubeMaterial()
{
    //same as glass material, no separate implementation needed for now.
}
void WOMMaterial::GenerateTubeInsideMaterial()
{
    G4NistManager* nist = G4NistManager::Instance();
    air = nist -> FindOrBuildMaterial("G4_AIR");

    std::vector<G4double> rIndexEnergy = { 1.775 * eV, 6.212 * eV};
    std::vector<G4double> rIndex = {1.0, 1.0};
    G4double numEntries = rIndexEnergy.size();

    mptAir = new G4MaterialPropertiesTable();
    mptAir -> AddProperty("RINDEX", rIndexEnergy, rIndex, numEntries);

    air -> SetMaterialPropertiesTable(mptAir);
}
std::vector<G4double>& WOMMaterial::WavelengthToEnergy(std::vector<G4double>& wLen)
{
    energy.clear();
    if(!wLen.empty())
    {
        for(int i = wLen.size() - 1; i >= 0; i--)
        {
            energy.push_back( (6.626 * pow(10, -34) * 3 * pow(10, 8)) / (wLen.at(i) * pow(10, -9) * 1.6 * pow(10, -19)));
        }
    }
    else
    {
        std::cerr << "Invalid array in WOMMaterial::WavelengthToEnergy(). Aborting..." << std::endl;
        abort();
    }

    return energy;
}
