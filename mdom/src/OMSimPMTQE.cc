///////////////////////////////////////////////////////////////////////////
//* ---------------------
//* J4PMTQEFunc Class
//* ---------------------
//* (Description)
//*   datum class for PMT Qe.
//*   Reads ascii quantum-efficiency data and interpolate between
//*   data points. Draw function provides TH1D histogram of QE vs WaveLen.
//*
//* (Requires)
//*     Ascii quantum-efficiency vs wavelength data
//*     See ReadQeTable() to refer proper format
//* (Provides)
//*	J4PMTQEFunc              : Qe modules
//* (Update Record)
//*	2004/8/14  S.Yoshida	Original version.
///////////////////////////////////////////////////////////////////////////

#include "OMSimPMTQE.hh"
#include "Interpolation.hh"
#include "G4SystemOfUnits.hh"

int kDegOfPolynominal = 6;
double kWAVELENGTH_MIN = (270. * nm);
double kWAVELENGTH_MAX = (720. * nm);
//Interpolation::fDelta = 1.0;
using namespace std;

//_____________________________________________________________________
OMSimPMTQE::OMSimPMTQE()
{
   fScale = 1.0;
}

//_____________________________________________________________________
OMSimPMTQE::~OMSimPMTQE()
{
}

//_____________________________________________________________________
/*void OMSimPMTQE::SetInputFileName(const G4String &fname)
{
   fQeDataName = fname;
}

//_____________________________________________________________________
void OMSimPMTQE::SetPMTModel(const G4String &pmtmodel)
{
   fQeDataName = "/home/waly/bulkice_doumeki/mdom/InputFile/TA0001_HamamatsuQE.data";
   //fQeDataName += pmtmodel;
   //fQeDataName += "/2D_QExCE/";
   //fQeDataName += pmtmodel;
   //fQeDataName += "_HamamatsuQE.data";
}
*/
//_____________________________________________________________________
void OMSimPMTQE::ReadQeTable()
{
//=====================================================================
//* Read the Qe data table ------------------------------------------
// data format must be :
//   lambda[m]  Qe[%]
// AND number of data is fixed to kNumber_of_data.
//=====================================================================

   Interpolation::DataReader(fQeDataName,
                             fWaveLengthTable, fQeTable);
   // original wavelength table is written in [nm].
   // Must scale to in geant4 default unit.
   int ndata = fWaveLengthTable.size();
    //std::cout << "++++++++++++++Size is " << ndata << " +++++++++" << std::endl;
   for (int i=0; i<ndata; i++) {
      fWaveLengthTable[i] *= nm;


      /*std::cerr << "wl : " << fWaveLengthTable[i]/nm
                << " qe : " << fQeTable[i]/nm << std::endl;
                */

   }


}

//_____________________________________________________________________
G4double OMSimPMTQE::GetQe(G4double lambda)
{
   // lambda must be in G4 unit (mm).
   // returns qe in percent.
   return Eval(lambda);
}

//_____________________________________________________________________
G4double OMSimPMTQE::Eval(G4double x)
{
//=====================================================================
// Get Quantum Efficieny
// x: Wave length in G4 unit [mm]
// The returned value is interpolated between data points
// using Interpolation class, then multiplied fScale.
// returns QE [%]
//=====================================================================
  double lambda = x;

  G4double qe= 0.0;
  int ndata = fWaveLengthTable.size();

  int mdeg;
  if(lambda < 290.0*nm) mdeg = 2;
  else mdeg = kDegOfPolynominal;
  qe = Interpolation::MThPolynominalInterpolate(fWaveLengthTable,
                                        fQeTable,
                                        ndata, lambda, mdeg);

  if (lambda<kWAVELENGTH_MIN || lambda> kWAVELENGTH_MAX) qe=0.0;
  if(qe<0.0) qe=0.0;

  return qe*fScale;

}
double OMSimPMTQE::RandomGen()
{
	std::random_device rd;  // Obtain a random seed from the hardware
    std::mt19937 gen(rd()); // Seed the generator

    std::uniform_real_distribution<double> dis(0.0, 1.0); // Define the distribution

    double random_num = dis(gen); // Generate a random number

    //std::cout << "Random number between 0 and 1: " << randomNum << std::endl;

    return random_num;
}
