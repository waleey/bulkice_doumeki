 //$Id: J4PMTQEFunc.h,v 1.3 2006/11/05 16:34:47 hoshina Exp $
#ifndef __OMSIMPMTQE__H__
#define __OMSIMPMTQE__H__

//*************************************************************************
//* ---------------------
//* J4PMTQEFunc Class
//* ---------------------
//* (Description)
//*     module class for PMT Qe
//* (Requires)
//* (Provides)
//*	J4PMTQEFunc              : Qe modules
//* (Update Record)
//*	2004/8/14  S.Yoshida	Original version.
//*************************************************************************
#include <iostream>
#include <fstream>
#include <vector>
//#include "J4VFunction.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include <random>

class OMSimPMTQE {

public:

  OMSimPMTQE();

   ~OMSimPMTQE();

  G4double GetQe(G4double waveLength);
  G4double GetScale()  const     { return fScale;}
  void             SetScale(G4double scale){ fScale = scale;}
  void             SetInputFileName(const G4String &fname);
  void             SetPMTModel(const G4String &pmtmodel);
  G4double Eval(G4double x);

  void ReadQeTable( );
  double RandomGen();

private:


  G4String  fQeDataName = "/home/waly/bulkice_doumeki/mdom/InputFile/TA0001_HamamatsuQE.data"; // name of qe data file

  std::vector<double> fQeTable;          // Qe data
  std::vector<double> fWaveLengthTable;  // Wavelength data

  // Renomalization Scale.
  // For accounting PMT individual variances.
  G4double fScale;                // scale factor

};

#endif
