#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "OMSimDataFileTypes.hh"
#include "OMSimInputData.hh"
#include "OMSimLogger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include <G4UnitsTable.hh>
#include <dirent.h>
#include <cmath>
#include <algorithm>
#include <vector>

namespace pt = boost::property_tree;
extern std::vector<double> readColumnDouble (G4String fn, int col);
void change_in_ascending_order(const std::vector<double> &ene, const std::vector<double> &val, std::vector<double> &ene2, std::vector<double> &val2)
{
    unsigned int n = ene.size();
    ene2.clear();
    ene2.resize(n);
    val2.clear();
    val2.resize(n);
    //std::cerr <<"vector size is " << n << std::endl;

    if (ene[0] < ene[n-1]) {
        //std::cerr << "do not need to make it reverse order" << std::endl;
        for (unsigned int i=0; i<n; ++i) {
            ene2[i] = ene[i];
            val2[i] = val[i];
            //std::cerr <<"photon energy, val = " << i << ", " << ene2[i] << ", " << val2[i] << std::endl;
        }
    } else {
        //std::cerr << "change to reverse order" << std::endl;
        for (unsigned int i=0; i<n; ++i) {
            ene2[i] = ene[n-1-i];
            val2[i] = val[n-1-i];
            //std::cerr <<"photon energy, val = " << i << ", " << ene2[i] << ", " << val2[i] << std::endl;
        }
    }
}


/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                                  Base Abstract Classes
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
/**
 * Constructor. Saves file name to member.
 * @param pFileName
 */
abcDataFile::abcDataFile(G4String pFileName)
{   mFileData = new ParameterTable();
    mFileName = pFileName;
}

template <typename T>
/**
 * Transforms the values inside a pTree-array to a vector. The values can be also transformed to a G4double.
 * @param pVector  vector where the (transformed) values are saved
 * @param pTree    pTree containing json data
 * @param pKey     json attribute label where values are found
 * @param pScaling Values of array are multiplied by this factor. You can set the fisical unit with this.
 * @param pInverse In case you need the inverse of a value x, 1/x is appended (e.g. transforming from nm to eV)
 */
void abcDataFile::ParseToVector(std::vector<T> &pVector, pt::ptree pTree, std::basic_string<char> pKey, G4double pScaling, bool pInverse)
{
    for (pt::ptree::value_type &ridx : pTree.get_child(pKey))
    { //get array from element with key "pKey" of the json
        if (pInverse)
        { // if we need 1/x
            pVector.push_back(pScaling / ridx.second.get_value<T>());
        }
        else
        { // otherwise we only by scaling factor
            pVector.push_back(ridx.second.get_value<T>() * pScaling);
        }
    }
}

/**
 * Defines new material from data in json-file.
 */
void abcMaterialData::CreateMaterial()
{
    //std::cerr <<"abcMaterialData::CreateMaterial is called" << std::endl;

    pt::read_json(mFileName, mJsonTree); //read json file into mJsonTree

    mFileData->AppendParameterTable(mFileName);
    mMPT = new G4MaterialPropertiesTable();
    mMatDatBase = G4NistManager::Instance();

    mObjectName = mJsonTree.get<G4String>("jName");
    const G4String lDataType = mJsonTree.get<G4String>("jDataType");

    const G4double lDensity = mFileData->GetValue(mObjectName, "jDensity");

    const G4String lState_str = mJsonTree.get<G4String>("jState");
    const G4State lState = GetState(lState_str);

    //Defining the material with its density, number of components, state and name
    mMaterial = new G4Material(mObjectName, lDensity, mJsonTree.get_child("jComponents").size(), lState);

    //Construct material with fractional components (isotopes or G4-Materials)
    for (pt::ptree::value_type &key : mJsonTree.get_child("jComponents"))
    {
        std::string componentName = key.first;
        double componentFraction = key.second.get_value<double>();
        mMaterial->AddMaterial(mMatDatBase->FindOrBuildMaterial(componentName), componentFraction);
    }
    G4String mssg = "New Material defined: " + mMaterial->GetName();
    info(mssg);

    //std::cerr <<"abcMaterialData::CreateMaterial end" << std::endl;
}
/**
 * Extracts absorption length and adds it to the material property table
 */
void abcMaterialData::ExtractAbsorptionLength()
{
    //std::cerr <<"abcMaterialData::ExtractAbsorptionLength is called" << std::endl;
    std::vector<G4double> lAbsLength, lAbsLength2;
    std::vector<G4double> lAbsLengthEnergy, lAbsLengthEnergy2;
    ParseToVector(lAbsLength, mJsonTree, "jAbsLength", 1 * mm, false);
    ParseToVector(lAbsLengthEnergy, mJsonTree, "jAbsLengthWavelength", mHC_eVnm, true);
    // must be in ascending order
    change_in_ascending_order(lAbsLengthEnergy, lAbsLength, lAbsLengthEnergy2, lAbsLength2);
    mMPT->AddProperty("ABSLENGTH", &lAbsLengthEnergy2[0], &lAbsLength2[0], static_cast<int>(lAbsLength.size()));
    //std::cerr <<"abcMaterialData::ExtractAbsorptionLength end" << std::endl;
}
/**
 * Extracts refraction index and adds it to the material property table
 */
 void abcMaterialData::ReverseCopy(std::vector<G4double>& arr1, G4double* arr2)
 {
    G4int len = arr1.size();
    for(int i = 0; i < len; i++)
    {
        arr2[len -1 - i] = arr1.at(i);
    }

 }
void abcMaterialData::ExtractRefractionIndex()
{
    //std::cerr <<"abcMaterialData::ExtractRefractionIndex is called" << std::endl;
    std::vector<G4double> lRefractionIndex, lRefractionIndex2;
    std::vector<G4double> lRefractionIndexEnergy, lRefractionIndexEnergy2;
    ParseToVector(lRefractionIndex, mJsonTree, "jRefractiveIdx", 1., false);
    ParseToVector(lRefractionIndexEnergy, mJsonTree, "jRefractiveIdxWavelength", mHC_eVnm, true);
    // must be in ascending order
    change_in_ascending_order(lRefractionIndexEnergy, lRefractionIndex, lRefractionIndexEnergy2, lRefractionIndex2);
    mMPT->AddProperty("RINDEX", &lRefractionIndexEnergy2[0], &lRefractionIndex2[0], static_cast<int>(lRefractionIndex.size()));
    //std::cerr <<"abcMaterialData::ExtractRefractionIndex end" << std::endl;
}
/**
*Setting up Scintillation Property for Vetrovex glass sample
*Vetrovex sample is referred as vas
*Spectrum taken from M.Unland's thesis
*input file containing the scintillation spectrum data is saved in the InputFile dir.
**/
void abcMaterialData::AddScintillationSpectrum(std::vector<G4double>& vasEnergy, std::vector<G4double>& vasScint)
{
    G4String filename = "../InputFile/VAS_Scintillation_Spectrum.data"; //Change it according to your file path
    //G4String filename = "/home/waly/bulkice_doumeki/mdom/InputFile/Vitrovex_scint.txt";
    //G4String filename = "/home/waly/bulkice_doumeki/mdom/InputFile/vis_corrected.txt";
    std::ifstream file(filename);
    if(!(file.is_open()))
    {
        std::cerr << "Error while opening " << filename << ". Aborting..." << std::endl;
        abort();
    }
    else
    {
        G4double en = 0;
        G4double intensity = 0;

        vasEnergy.clear();
        vasScint.clear();

        while(file >> en >> intensity)
        {
            vasEnergy.push_back(en * eV);
            vasScint.push_back(intensity);
        }
        file.close();
        file.clear();
    }

}
/**
*Adding Scintillation property to materials. It will only be applied to the Pressure Vessel Glass Material.
*For now, dummy scintillation yield, and scintillation decay time constants are used.
*Particle independent scintillation is assumed for starting. Will change it later.
*Scintillation time constant is the average scintillation time constant
*for Vetrovex glass after correction at -20C temperature.
*Data taken from M. Dittmer's MA thesis.
*Alpha Scintillation Yield is 75.0
*Electron yield is 129 (2018 batch sample)
*Scintillation is only defined for alpha and electron.
*Yield for other particles are unknown.
**/
/*void abcMaterialData::AddScintillationProperty()
{
    std::vector<G4double> vasEnergy;
    std::vector<G4double> vasScint;
    AddScintillationSpectrum(vasEnergy, vasScint);

    std::vector<G4double> alphaEnergy = {0, 100 * MeV};
    std::vector<G4double> alphaYield = {5, 75 * 100};
    std::vector<G4double> electronEnergy = {0, 100 * MeV};
    std::vector<G4double> electronYield = {5, 129 * 100};
    G4double generalYield = 129; //set the yield for electron

    mMPT->AddProperty("SCINTILLATIONCOMPONENT1", vasEnergy, vasScint);
    mMPT->AddConstProperty("SCINTILLATIONYIELD", generalYield/MeV);
    mMPT->AddProperty("ALPHASCINTILLATIONYIELD", alphaEnergy, alphaYield);
    mMPT->AddConstProperty("ALPHASCINTILLATIONYIELD1", 1.0);
    mMPT->AddProperty("ELECTRONSCINTILLATIONYIELD" , electronEnergy, electronYield);
    //mMPT->AddProperty("POSITRONSCINTILLATIONYIELD", electronEnergy, electronYield);
    mMPT->AddConstProperty("ELECTRONSCINTILLATIONYIELD1", 1.0);
    mMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 8.83 * microsecond);

    //mMaterial->GetIonisation()->SetBirksConstant(0.0151*cm/MeV); commenting out Birks constant cuz it's not defined in any previous lit.
}*/
/**
 * State in string to G4State
 * @param  G4String
 * @return G4State
 */
G4State abcMaterialData::GetState(G4String pState_str)
{
    G4State lState;
    if (pState_str == "kStateSolid")
        lState = kStateSolid;
    else if (pState_str == "kStateLiquid")
        lState = kStateLiquid;
    else if (pState_str == "kStateGas")
        lState = kStateGas;
    else
        lState = kStateUndefined;
    return lState;
}
/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                                     Derived Classes
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
/**
 * Extracts and creates material for material with refraction index and absorption length defined.
 */
void RefractionAndAbsorption::ExtractInformation()
{
    CreateMaterial();
    ExtractAbsorptionLength();
    ExtractRefractionIndex();
    /*if(fisGlass)
    {
        std::cout << "********Scintillation Property Activated*******" << std::endl;
        AddScintillationProperty();
    }*/
    if (mMaterial -> GetName() == "RiAbs_Glass_Tube"){
        //critical("Found");
        std::cout << "Adding Sctintillation for Vitro Glass Tube..." << std::endl;
        G4double factor = 0.2;
        G4double times20[7] = {
        1.500932040523042781e-03*s*factor*1.5,
        4e-02*s,
        2.304116799105801212e-04*s*factor*1.5,
        2.967754892640673745e-07*s*factor,
        2.193657986806030323e-06*s*factor,
        4.666464222616022482e-05*s*factor*1.5,
        1.012137273486919956e-05*s*factor
        /**
            4e-02*s,
        1.500932040523042781e-03*s*factor*1.5,
        2.304116799105801212e-04*s*factor*1.5,
        4.666464222616022482e-05*s*factor*1.5,
        1.012137273486919956e-05*s*factor,
        2.193657986806030323e-06*s*factor,
        2.967754892640673745e-07*s*factor*/
    };
    G4double amplitudes20[7] = {

    6.290979060619150687e-02*0.6,
    4e-02,
    6.846028776435156282e-02,
    1.630113195117202096e-01*0.5,
    3.496767714671291660e-01*0.35,
    1.326011725433764721e-01,
    2.233406581072310548e-01
    /**
        4e-02,
        6.290979060619150687e-02*0.6,
        6.846028776435156282e-02,
        1.326011725433764721e-01,
        2.233406581072310548e-01,
        3.496767714671291660e-01*0.35,
        1.630113195117202096e-01*0.5*/
    };
    //std::sort(amplitudes20, amplitudes20 + sizeof(amplitudes20) / sizeof(amplitudes20[0]));
    //G4String DataFile = "../InputFile/VAS_Scintillation_Spectrum.data";
    G4String DataFile = "../InputFile/Vitrovex_scint.txt"; //will change soon
    std::vector<double> fileFirstColumn = readColumnDouble(DataFile, 1);
   /* for(auto value : fileFirstColumn)
    {
        std::cout << value << std::endl;
    }*/
    std::vector<double> fileSecondColumn = readColumnDouble(DataFile, 2);

    G4double VV_WL[int(fileFirstColumn.size())];
    G4double VV_I[int(fileSecondColumn.size())];
    //std::copy_backward(fileFirstColumn.begin(), fileFirstColumn.end(),  VV_WL.end());
    //std::copy_backward(fileSecondColumn.begin(), fileSecondColumn.end(), VV_I.end());
    ReverseCopy(fileFirstColumn, VV_WL);
    ReverseCopy(fileSecondColumn, VV_I);
    for (unsigned int u = 0; u <fileFirstColumn.size(); u++) {

        VV_WL[u] = 1239.84193 / VV_WL[u]*eV;

    }


    mMPT->AddProperty("FIRSTCOMPONENT", VV_WL,VV_I,138, true, true);
    mMPT->AddConstProperty("SCINTILLATIONYIELD",59.6/MeV);
    mMPT->AddProperty("FRACTIONLIFETIMES",amplitudes20,times20,7, true, true);
    mMPT->AddProperty("SCINTILLATIONSPECTRUM", VV_WL,VV_I,138, true, true);
    mMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    }


    if (mMaterial->GetName() == "RiAbs_Glass_Vitrovex"){
        //critical("Found");
        std::cout << "Adding Sctintillation for Vitro nGlass Vessel..." << std::endl;
       // G4double factor = 0.3;
        G4double times20[6] = {
        1.500932040523042781e-03*s,
        2.304116799105801212e-04*s,
        4.666464222616022482e-05*s,
        2.967754892640673745e-07*s,
        1.012137273486919956e-05*s,
        2.193657986806030323e-06*s
    };
    G4double amplitudes20[6] = {
        6.290979060619150687e-02,
        6.846028776435156282e-02,
        1.326011725433764721e-01,
        1.630113195117202096e-01,
        2.233406581072310548e-01,
        3.496767714671291660e-01
    };
    //G4String DataFile = "../InputFile/VAS_Scintillation_Spectrum.data";
    G4String DataFile = "../InputFile/Vitrovex_scint.txt";
    std::vector<double> fileFirstColumn = readColumnDouble(DataFile, 1);
    std::vector<double> fileSecondColumn = readColumnDouble(DataFile, 2);

    G4double VV_WL[int(fileFirstColumn.size())];
    G4double VV_I[int(fileSecondColumn.size())];
    ReverseCopy(fileFirstColumn, VV_WL);
    ReverseCopy(fileSecondColumn, VV_I);
    //std::copy_backward(fileFirstColumn.begin(), fileFirstColumn.end(), std::end(VV_WL));
    //std::copy_backward(fileSecondColumn.begin(), fileSecondColumn.end(), std::end(VV_I));
    /*for(auto value : VV_WL)
    {
        std::cout << value << std::endl;
    }*/
    for (unsigned int u = 0; u <fileFirstColumn.size(); u++) {

        VV_WL[u] = 1239.84193 / VV_WL[u]*eV;

    }

   mMPT->AddProperty("FIRSTCOMPONENT",VV_WL,VV_I,138, true, true);
    mMPT->AddConstProperty("SCINTILLATIONYIELD",50/MeV);
    mMPT->AddProperty("FRACTIONLIFETIMES",amplitudes20,times20,6, true, true);
    mMPT->AddProperty("SCINTILLATIONSPECTRUM",VV_WL,VV_I,138, true, true);
    mMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    }


   /* if (mMaterial->GetName() == "RiAbs_Glass_Chiba"){
        //critical("Found");
        std::cout << "Adding Sctintillation for Chiba Glass..." << std::endl;
        G4double DEGGtimes20[5] = {
3.7596973315004064e-05*s,
8.480364129990657e-06*s,
2.0193372706901436e-06*s,
1.7722047480973938e-07*s,
0.00024852344312136036*s
    };
        G4double DEGGamplitudes20[5] = {
0.2644264442708002,
0.15726554891967992,
0.209532534040833,
0.1780995886502345,
0.19067588411845232
    };
    G4String DataFile = "../InputFile/Vitrovex_scint.txt";
    std::vector<double> fileFirstColumn = readColumnDouble(DataFile, 1);
    std::vector<double> fileSecondColumn = readColumnDouble(DataFile, 2);

    G4double VV_WL[138];
    G4double VV_I[138];
    ReverseCopy(fileFirstColumn, VV_WL);
    ReverseCopy(fileSecondColumn, VV_I);
    //std::copy_backward(fileFirstColumn.begin(), fileFirstColumn.end(), std::end(VV_WL));
    //std::copy_backward(fileSecondColumn.begin(), fileSecondColumn.end(), std::end(VV_I));

    for (unsigned int u = 0; u <fileFirstColumn.size(); u++) {

        VV_WL[u] = 1239.84193 / VV_WL[u]*eV;

    }

    mMPT->AddProperty("FIRSTCOMPONENT",VV_WL,VV_I,138);
    mMPT->AddConstProperty("SCINTILLATIONYIELD",55/MeV);
    mMPT->AddProperty("FRACTIONLIFETIMES",DEGGamplitudes20,DEGGtimes20,5);
    mMPT->AddProperty("SCINTILLATIONSPECTRUM",VV_WL,VV_I,138);
    mMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    }*/
    mMaterial->SetMaterialPropertiesTable(mMPT);
}
/**
 * Extracts and creates material for material with refraction index defined.
 */
void RefractionOnly::ExtractInformation()
{
    CreateMaterial();
    ExtractRefractionIndex();
    mMaterial->SetMaterialPropertiesTable(mMPT);
}
/**
 * Extracts and creates material without optical properties.
 */
void NoOptics::ExtractInformation()
{
    CreateMaterial();
}
/**
 * Extracts and creates ice with optical properties from IceCube.
 */
void IceCubeIce::ExtractInformation()
{
    CreateMaterial(); //creates IceCubeICE

    G4Material *lIceMie = new G4Material("IceCubeICE_SPICE", mFileData->GetValue(mObjectName,"jDensity"), mMatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid); //create IceCubeICE_SPICE
    G4Material *lBubleColumnMie = new G4Material("Mat_BubColumn", mFileData->GetValue(mObjectName,"jDensity"), mMatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid); //create IceCubeICE_SPICE
    std::vector<G4double> lMieScatteringLength, lMieScatteringLength2;
    std::vector<G4double> lMieScatteringLength_BubleColumn, lMieScatteringLength_BubleColumn2;
    std::vector<G4double> lWavelength;
    std::vector<G4double> lRefractionIndex, lRefractionIndex2;
    std::vector<G4double> lRefractionIndexEnergy, lRefractionIndexEnergy2;
    std::vector<G4double> lAbsLength, lAbsLength2;
    ParseToVector(lRefractionIndexEnergy, mJsonTree, "jWavelength_spice", mHC_eVnm, true);
    ParseToVector(lWavelength, mJsonTree, "jWavelength_spice", 1 * nm, false);
    ParseToVector(mSpice_be400inv, mJsonTree, "jbe400inv_spice", 1 * m, false);
    ParseToVector(mSpice_a400inv, mJsonTree, "ja400inv_spice", 1 * m, false);
    ParseToVector(mSpice_Depth, mJsonTree, "jDepth_spice", 1 * m, false);

    for (int u = 0; u < static_cast<int>(lRefractionIndexEnergy.size()); u++)
    {
        lRefractionIndex.push_back(Spice_Refraction(lWavelength.at(u)));
        lAbsLength.push_back(Spice_Absorption(lWavelength.at(u)));
        lMieScatteringLength.push_back(Mie_Scattering(lWavelength.at(u)));
        lMieScatteringLength_BubleColumn.push_back(mInnercolumn_b_inv);
    }

    // new G4 requires energies in ascending order.
    change_in_ascending_order(lRefractionIndexEnergy, lRefractionIndex, lRefractionIndexEnergy2, lRefractionIndex2);
    change_in_ascending_order(lRefractionIndexEnergy, lAbsLength, lRefractionIndexEnergy2, lAbsLength2);
    change_in_ascending_order(lRefractionIndexEnergy, lMieScatteringLength, lRefractionIndexEnergy2, lMieScatteringLength2);
    change_in_ascending_order(lRefractionIndexEnergy, lMieScatteringLength_BubleColumn, lRefractionIndexEnergy2, lMieScatteringLength_BubleColumn2);

    //give refractive index to IceCubeICE. This is used also for IceCubeICE_SPICE
    mMPT->AddProperty("RINDEX", &lRefractionIndexEnergy2[0], &lRefractionIndex2[0], static_cast<int>(lRefractionIndex.size()));
    mMaterial->SetMaterialPropertiesTable(mMPT);

    //give properties to IceCubeICE_SPICE
    G4MaterialPropertiesTable* lMPT_spice = new G4MaterialPropertiesTable();
    lMPT_spice->AddProperty("RINDEX", &lRefractionIndexEnergy2[0], &lRefractionIndex2[0], static_cast<int>(lRefractionIndex.size()));
    lMPT_spice->AddProperty("ABSLENGTH", &lRefractionIndexEnergy2[0], &lAbsLength2[0], static_cast<int>(lAbsLength.size()));
    lMPT_spice->AddProperty("MIEHG", &lRefractionIndexEnergy2[0], &lMieScatteringLength2[0], static_cast<int>(lRefractionIndex.size()));
    lMPT_spice->AddConstProperty("MIEHG_FORWARD", mMIE_spice_const[0]);
    lMPT_spice->AddConstProperty("MIEHG_BACKWARD", mMIE_spice_const[1]);
    lMPT_spice->AddConstProperty("MIEHG_FORWARD_RATIO", mMIE_spice_const[2]);
    lIceMie->SetMaterialPropertiesTable(lMPT_spice);
    G4String mssg = "Ice properties at depth " + std::to_string(mSpice_Depth[mSpiceDepth_pos] / m) + " m.";
    notice(mssg);
    //now give the properties to the bubble column, which are basically the same ones but with the chosen scattering lenght
    G4MaterialPropertiesTable* lMPT_holeice = new G4MaterialPropertiesTable();
    lMPT_holeice->AddProperty("RINDEX", &lRefractionIndexEnergy2[0], &lRefractionIndex2[0], static_cast<int>(lRefractionIndex.size()));
    lMPT_holeice->AddProperty("ABSLENGTH", &lRefractionIndexEnergy2[0], &lAbsLength2[0], static_cast<int>(lAbsLength.size()));
    lMPT_holeice->AddProperty("MIEHG", &lRefractionIndexEnergy2[0], &lMieScatteringLength_BubleColumn2[0], static_cast<int>(lRefractionIndex.size()));
    lMPT_holeice->AddConstProperty("MIEHG_FORWARD", mMIE_spice_const[0]);
    lMPT_holeice->AddConstProperty("MIEHG_BACKWARD", mMIE_spice_const[1]);
    lMPT_holeice->AddConstProperty("MIEHG_FORWARD_RATIO", mMIE_spice_const[2]);
    lBubleColumnMie->SetMaterialPropertiesTable(lMPT_holeice);
}
/*
 * %%%%%%%%%%%%%%%% Functions for icecube ice optical properties %%%%%%%%%%%%%%%%
 */
/**
 * This gives you temperature of ice depending on the depth.
 * Function needed for the calculation of scattering and absorption length of the ice.
 * @param pDepth Depth in m from where we need the temperature
 * @return Temperature
 */
G4double IceCubeIce::Spice_Temperature(G4double pDepth)
{
    G4double spice_temp = 221.5 - 0.00045319 / m * pDepth + 5.822e-6 / m2 * pow(pDepth, 2.);
    return spice_temp;
}

/**
 * Calculation of the absorption length of IC-ice for a specific wavelength
 * @param pLambd Wavelength
 * @return Absorption length
 */
G4double IceCubeIce::Spice_Absorption(G4double pLambd)
{
    G4double lKappa = 1.08;
    G4double lParamA = 6954. / m;
    G4double lParamB = 6618 * nm;
    G4double lAdust = 1. / (mSpice_a400inv[mSpiceDepth_pos]) * pow(pLambd / (400. * nm), -lKappa);
    G4double lDeltaTau = Spice_Temperature(mSpice_Depth[mSpiceDepth_pos]) - Spice_Temperature(1730.);
    G4double la_inv = 1. / (lAdust + lParamA * exp(-lParamB / pLambd) * (1. + 0.01 * lDeltaTau));
    return la_inv;
}

/**
 * Calculation of the refraction index of IC-ice for a specific wavelength.
 * @param pLambd Wavelength
 * @return Refraction index
 */
G4double IceCubeIce::Spice_Refraction(G4double pLambd)
{
    // unknown depth. Parametrization by Thomas Kittler.
    G4double lLambd3 = pLambd * 1e-3;
    G4double lNphase = 1.55749 - 1.57988 / nm * lLambd3 + 3.99993 / (nm * nm) * pow(lLambd3, 2) - 4.68271 / (nm * nm * nm) * pow(lLambd3, 3) + 2.09354 / (nm * nm * nm * nm) * pow(lLambd3, 4);
    return lNphase; // using this now after discussion with Timo
}
/**
 * Calculation of the mie scattering length of IC-ice for a specific wavelength
 * @param pLambd Wavelength
 * @return Mie scattering length
 */
G4double IceCubeIce::Mie_Scattering(G4double pLambd)
{
    // depth_pos is the coordinate for the chosen depth in Depth_spice. For example to choose
    // depth=2278.2 m, we use depth_pos = 88
    G4double lAlpha = 0.90;
    G4double lAv_costheta = 0.9;
    G4double lBe_inv = 1. / (1. / (mSpice_be400inv[mSpiceDepth_pos]) * pow((pLambd / (400. * nm)), -lAlpha));
    G4double lB_inv = lBe_inv * (1. - lAv_costheta);
    return lB_inv;
}
/*
 * %%%%%%%%%%%%%%%% Functions of derived class ReflectiveSurface %%%%%%%%%%%%%%%%
 */
/**
 * Defines new reflective surface from data in json-file.
 */
void ReflectiveSurface::ExtractInformation()
{

    pt::read_json(mFileName, mJsonTree); //read json file into mJsonTree

    mObjectName = mJsonTree.get<G4String>("jName");
    G4String lModelStr = mJsonTree.get<G4String>("jModel");
    G4String lFinishStr = mJsonTree.get<G4String>("jFinish");
    G4String lTypeStr = mJsonTree.get<G4String>("jType");

    G4OpticalSurfaceModel lModel = GetOpticalSurfaceModel(lModelStr);
    G4OpticalSurfaceFinish lFinish = GetOpticalSurfaceFinish(lFinishStr);
    G4SurfaceType lType = GetSurfaceType(lTypeStr);

    mOpticalSurface = new G4OpticalSurface(mObjectName, lModel, lFinish, lType);
    G4MaterialPropertiesTable *lMPT = new G4MaterialPropertiesTable();

    try // Only few materials have jSigmaAlpha defined
    {
        G4double lSigmaAlpha = mJsonTree.get<G4double>("jSigmaAlpha");
        mOpticalSurface->SetSigmaAlpha(lSigmaAlpha);
    }
    catch (...)
    {
    } // not very elegant, I know...

    for (pt::ptree::value_type &key : mJsonTree.get_child("jProperties"))
    {
        G4String lKey = key.second.get_value<G4String>();
        std::vector<G4double> lPhotonEnergy, lPhotonEnergy2;
        std::vector<G4double> lValues, lValues2;
        ParseToVector(lValues, mJsonTree, "jValues_" + lKey, 1., false);
        ParseToVector(lPhotonEnergy, mJsonTree, "jWavelength_" + lKey, mHC_eVnm, true);
        // energy must be in ascending order
        change_in_ascending_order(lPhotonEnergy, lValues, lPhotonEnergy2, lValues2);
        lMPT->AddProperty(lKey, &lPhotonEnergy2[0], &lValues2[0], static_cast<int>(lPhotonEnergy.size()));
    }

    mOpticalSurface->SetMaterialPropertiesTable(lMPT);
    G4String mssg = "New Optical Surface: " + mObjectName;
    info(mssg);
}

/**
 * OpticalSurfaceFinish in string to G4OpticalSurfaceFinish
 * @param  G4String
 * @return G4OpticalSurfaceFinish
 */
G4OpticalSurfaceFinish ReflectiveSurface::GetOpticalSurfaceFinish(G4String pFinish)
{
    G4OpticalSurfaceFinish lFinish;
    if (pFinish == "polished")
        lFinish = polished;
    else if (pFinish == "polishedfrontpainted")
        lFinish = polishedfrontpainted;
    else if (pFinish == "polishedbackpainted")
        lFinish = polishedbackpainted;
    else if (pFinish == "ground")
        lFinish = ground;
    else if (pFinish == "groundfrontpainted")
        lFinish = groundfrontpainted;
    else if (pFinish == "groundbackpainted")
        lFinish = groundbackpainted;
    else if (pFinish == "polishedlumirrorair")
        lFinish = polishedlumirrorair;
    else if (pFinish == "polishedlumirrorglue")
        lFinish = polishedlumirrorglue;
    else if (pFinish == "polishedair")
        lFinish = polishedair;
    else if (pFinish == "polishedteflonair")
        lFinish = polishedteflonair;
    else if (pFinish == "polishedtioair")
        lFinish = polishedtioair;
    else if (pFinish == "polishedtyvekair")
        lFinish = polishedtyvekair;
    else if (pFinish == "polishedvm2000air")
        lFinish = polishedvm2000air;
    else if (pFinish == "polishedvm2000glue")
        lFinish = polishedvm2000glue;
    else if (pFinish == "etchedlumirrorair")
        lFinish = etchedlumirrorair;
    else if (pFinish == "etchedlumirrorglue")
        lFinish = etchedlumirrorglue;
    else if (pFinish == "etchedair")
        lFinish = etchedair;
    else if (pFinish == "etchedteflonair")
        lFinish = etchedteflonair;
    else if (pFinish == "etchedtioair")
        lFinish = etchedtioair;
    else if (pFinish == "etchedtyvekair")
        lFinish = etchedtyvekair;
    else if (pFinish == "etchedvm2000air")
        lFinish = etchedvm2000air;
    else if (pFinish == "etchedvm2000glue")
        lFinish = etchedvm2000glue;
    else if (pFinish == "groundlumirrorair")
        lFinish = groundlumirrorair;
    else if (pFinish == "groundlumirrorglue")
        lFinish = groundlumirrorglue;
    else if (pFinish == "groundair")
        lFinish = groundair;
    else if (pFinish == "groundteflonair")
        lFinish = groundteflonair;
    else if (pFinish == "groundtioair")
        lFinish = groundtioair;
    else if (pFinish == "groundtyvekair")
        lFinish = groundtyvekair;
    else if (pFinish == "groundvm2000air")
        lFinish = groundvm2000air;
    else if (pFinish == "groundvm2000glue")
        lFinish = groundvm2000glue;
    return lFinish;
}
/**
 * OpticalSurfaceModel in string to G4OpticalSurfaceModel
 * @param  G4String
 * @return G4OpticalSurfaceModel
 */
G4OpticalSurfaceModel ReflectiveSurface::GetOpticalSurfaceModel(G4String pModel)
{
    G4OpticalSurfaceModel lModel;
    if (pModel == "glisur")
        lModel = glisur;
    else if (pModel == "unified")
        lModel = unified;
    else if (pModel == "LUT")
        lModel = LUT;
    return lModel;
}
/**
 * SurfaceType in string to G4SurfaceType
 * @param  G4String
 * @return G4SurfaceType
 */
G4SurfaceType ReflectiveSurface::GetSurfaceType(G4String pType)
{

    G4SurfaceType lType;
    if (pType == "dielectric_metal")
        lType = dielectric_metal;
    else if (pType == "dielectric_dielectric")
        lType = dielectric_dielectric;
    else if (pType == "dielectric_LUT")
        lType = dielectric_LUT;
    else if (pType == "firsov")
        lType = firsov;
    else if (pType == "x_ray")
        lType = x_ray;
    return lType;
}

