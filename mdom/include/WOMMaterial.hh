#ifndef WOMMATERICAL_HH_INCLUDED
#define WOMMATERICAL_HH_INCLUDED
/**
*Very rough draft of possible WOM Simulation
*Waly M Z Karim
*8/30/2023
**/

#include "G4Types.hh"
#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"

class WOMMaterial
{
public:
    WOMMaterial();
    ~WOMMaterial();

    inline G4Material* GetGlassMaterial() { return quartz; }
    inline G4Material* GetFillerMaterial() { return filler; }
    inline G4Material* GetPaintMaterial() { return quartzPaint; }
    inline G4Material* GetTubeMaterial() { return quartz; }
    inline G4Material* GetTubeInsideMaterial() { return air; }

private:
    void ReadData(const G4String&, std::vector<G4double>&);
    void GetSharedData(); //should be called in GenerateMaterial()
    void GenerateMaterials();
    void GenerateGlassMaterial();
    void GenerateFillerMaterial();
    void GeneratePaintMaterial();
    void GenerateTubeMaterial();
    void GenerateTubeInsideMaterial();
    std::vector<G4double>& WavelengthToEnergy(std::vector<G4double>&);

    G4Material* quartz;
    G4Material* quartzPaint;
    G4Material* filler;
    G4Material* air;

    G4MaterialPropertiesTable* mptQuartz;
    G4MaterialPropertiesTable* mptQuartzPaint;
    G4MaterialPropertiesTable* mptFiller;
    G4MaterialPropertiesTable* mptAir;

    G4double paintThickness;
    G4double quartzDensity;
    G4double fillerDensity;
    std::vector<G4double> wLen;
    std::vector<G4double> energy;

    //for defining materials
    G4int z, nElements;
    G4double a;

};
#endif // WOMMATERICAL_HH_INCLUDED
