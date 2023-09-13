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

    std::vector<G4double> rIndexEnergy =
    {1.775 * eV, 1.8736 * eV, 1.9722 * eV, 2.0708 * eV, 2.1694 * eV, 2.268* eV,  2.3666* eV, 2.4652* eV,
       2.5638* eV, 2.6624* eV, 2.761* eV,  2.8596* eV, 2.9582* eV, 3.0568* eV, 3.1554* eV, 3.254* eV,
       3.3526* eV, 3.4512* eV, 3.5498* eV, 3.6484* eV, 3.747* eV , 3.8456* eV, 3.9442* eV, 4.0428* eV,
       4.1414* eV, 4.24* eV  , 4.3386* eV, 4.4372* eV, 4.5358* eV, 4.6344* eV, 4.733* eV , 4.8316* eV,
       4.9302* eV, 5.0288* eV, 5.1274* eV, 5.226* eV , 5.3246* eV, 5.4232* eV, 5.5218* eV, 5.6204* eV,
       5.719* eV , 5.8176* eV, 5.9162* eV, 6.0148* eV, 6.1134* eV, 6.212* eV };

    G4int numEntries = rIndexEnergy.size();
    std::vector<G4double> rIndex;
    std::vector<G4double> absLen; //dummy value
    for(int i = 0; i < numEntries; i++)
    {
        rIndex.push_back(1.46);
        absLen.push_back(1.0 * m);
    }

    mptQuartz = new G4MaterialPropertiesTable();
    mptQuartz -> AddProperty("RINDEX", rIndexEnergy, rIndex, numEntries);
    mptQuartz -> AddProperty("ABSLENGTH", rIndexEnergy, absLen, numEntries);

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

    std::vector<G4double> rIndexEnergy =
    {1.775 * eV, 1.8736 * eV, 1.9722 * eV, 2.0708 * eV, 2.1694 * eV, 2.268* eV , 2.3666* eV, 2.4652* eV,
       2.5638* eV, 2.6624* eV, 2.761* eV,  2.8596* eV, 2.9582* eV, 3.0568* eV, 3.1554* eV, 3.254* eV,
       3.3526* eV, 3.4512* eV, 3.5498* eV, 3.6484* eV, 3.747* eV , 3.8456* eV, 3.9442* eV, 4.0428* eV,
       4.1414* eV, 4.24* eV  , 4.3386* eV, 4.4372* eV, 4.5358* eV, 4.6344* eV, 4.733* eV , 4.8316* eV,
       4.9302* eV, 5.0288* eV, 5.1274* eV, 5.226* eV , 5.3246* eV, 5.4232* eV, 5.5218* eV, 5.6204* eV,
       5.719* eV , 5.8176* eV, 5.9162* eV, 6.0148* eV, 6.1134* eV, 6.212* eV };
    G4double numEntries = rIndexEnergy.size();
    std::vector<G4double> rIndex;
    for(int i = 0; i < numEntries; i++)
    {
        rIndex.push_back(1.33);
        //absLen.push_back(1.0 * m);
    }

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

    std::vector<G4double> rIndexEnergy =
    {1.775 * eV, 1.8736 * eV, 1.9722 * eV, 2.0708 * eV, 2.1694 * eV, 2.268* eV,  2.3666* eV, 2.4652* eV,
       2.5638* eV, 2.6624* eV, 2.761* eV,  2.8596* eV, 2.9582* eV, 3.0568* eV, 3.1554* eV, 3.254* eV,
       3.3526* eV, 3.4512* eV, 3.5498* eV, 3.6484* eV, 3.747* eV , 3.8456* eV, 3.9442* eV, 4.0428* eV,
       4.1414* eV, 4.24* eV  , 4.3386* eV, 4.4372* eV, 4.5358* eV, 4.6344* eV, 4.733* eV , 4.8316* eV,
       4.9302* eV, 5.0288* eV, 5.1274* eV, 5.226* eV , 5.3246* eV, 5.4232* eV, 5.5218* eV, 5.6204* eV,
       5.719* eV , 5.8176* eV, 5.9162* eV, 6.0148* eV, 6.1134* eV, 6.212* eV };
    std::vector<G4double> rIndex;
    for(int i = 0; i < rIndexEnergy.size(); i++)
    {
        rIndex.push_back(1.46);
    }

    std::vector<G4double> absIntensity;
    std::vector<G4double> absLength;
    std::vector<G4double> absEnergy;
    ReadData("../InputFile/Absorption_WLS.txt", absIntensity);

    absEnergy = WavelengthToEnergy(wLen);

    for(int i = absIntensity.size() -1; i >= 0; i--)
    {
        absLength.push_back(- paintThickness / (log(1 - absIntensity.at(i) / 100)));
        //rIndexEnergy.push_back(absEnergy.at(i));
        //rIndex.push_back(1.46);
    }

    //rIndexEnergy.push_back(7.0 * eV);
    //rIndex.push_back(1.46);

    for(int i = 0; i < absLength.size(); i++)
    {
        std::cout << "absLen: " << absLength.at(i) / micrometer << std::endl;
        std::cout << "Energy: " << absEnergy.at(i) / eV << std::endl;
    }
    //exit(0);
    absIntensity.clear();

    std::vector<G4double> emIntensityInverted;
    std::vector<G4double> emIntensity;
    std::vector<G4double> emEnergy;

    ReadData("../InputFile/Emission_WLS.txt", emIntensityInverted);
    for(int i = emIntensityInverted.size() -1; i >= 0; i--)
    {
        emIntensity.push_back(emIntensityInverted.at(i) / 100);
    }
    emEnergy = WavelengthToEnergy(wLen);
    for(int i = 0; i < emIntensity.size(); i++)
    {
        std::cout << "EM: " << emIntensity.at(i) << std::endl;
        std::cout << "energy: " << emEnergy.at(i) / eV << std::endl;
    }
    //exit(0);
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

    std::vector<G4double> rIndexEnergy =
     {1.775 * eV, 1.8736 * eV, 1.9722 * eV, 2.0708 * eV, 2.1694 * eV, 2.268* eV,  2.3666* eV, 2.4652* eV,
       2.5638* eV, 2.6624* eV, 2.761* eV,  2.8596* eV, 2.9582* eV, 3.0568* eV, 3.1554* eV, 3.254* eV,
       3.3526* eV, 3.4512* eV, 3.5498* eV, 3.6484* eV, 3.747* eV , 3.8456* eV, 3.9442* eV, 4.0428* eV,
       4.1414* eV, 4.24* eV  , 4.3386* eV, 4.4372* eV, 4.5358* eV, 4.6344* eV, 4.733* eV , 4.8316* eV,
       4.9302* eV, 5.0288* eV, 5.1274* eV, 5.226* eV , 5.3246* eV, 5.4232* eV, 5.5218* eV, 5.6204* eV,
       5.719* eV , 5.8176* eV, 5.9162* eV, 6.0148* eV, 6.1134* eV, 6.212* eV };
    std::vector<G4double> rIndex;
    G4double numEntries = rIndexEnergy.size();
    for(int i = 0; i < numEntries; i++)
    {
        rIndex.push_back(1.0);
    }

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
            energy.push_back(( (6.626 * pow(10, -34) * 3 * pow(10, 8)) / (wLen.at(i) * pow(10, -9) * 1.6 * pow(10, -19))) * eV);
        }
    }
    else
    {
        std::cerr << "Invalid array in WOMMaterial::WavelengthToEnergy(). Aborting..." << std::endl;
        abort();
    }

    return energy;
}
