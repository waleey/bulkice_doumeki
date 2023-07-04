#ifndef __Interpolation__H__
#define __Interpolation__H__

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

#include "G4String.hh"
#include "G4SystemOfUnits.hh"


class Interpolation {

  private:
  static double fDelta;

 public:
  static int    SearchIndex(const std::vector<double>& xx,
                            const double x, const int n)
    {
    //******************************************************************
    // Search the array index such that x is between xx[index]
    // and xx[index+1]. xx must be monotonic, either increasing
    // or decreasing. index =0 or index = n-1 is returned
    // to indicate x is out of range.  n gives number of the elements
    // i.e. xx[0....n-1].
    //******************************************************************
      int jm;
      int jl=0;
      int ju = n;
      int index;
      bool ascnd;

      ascnd=(xx[n-1] >= xx[0]);
      while (ju-jl > 1) {
        jm=(ju+jl) >> 1;
        if ((x >= xx[jm]) == ascnd)
          jl=jm;
        else
          ju=jm;
      }
      if (x == xx[0]) index=0;
      else if(x == xx[n-1]) index=n-1;
      else index=jl;

      return index;
    }

  static double PolynominalInterpolate(const std::vector<double>& xa,
					 const std::vector<double>& ya,
					 const int n,
					 const double x)
    {
    //******************************************************************
    // Given arrays xa[0....n-1] and ya[0....n-1], and given a value x,
    // this method returns a value y by the polynominal interpolation.
    // Number of data ponts n is qeual to the plynominal of degree.
    //******************************************************************
      int dim = xa.size();
      if(dim>n) dim = n;
      std::vector<double> c(ya);
      std::vector<double> d(ya);

      double dif = fabs(x-xa[0]);
      double difNew;
      int closestIndex = 0;

      for(int i=0;i<dim;i++){ // Here we find the index of the closest table.
        difNew = fabs(x-xa[i]);
        if(difNew< dif){
          closestIndex = i;
          dif = difNew;
        }
      }

      double yOut = ya[closestIndex--];
      for(int mi=1; mi<dim; mi++){
        for(int i=0; i<dim - mi; i++){
          double ho = xa[i]-x;
          double hp = xa[i + mi]-x;
          double w = c[i+1]- d[i];
          double den = ho - hp;
          if(den == 0.0){
            std::cerr << "Error in polynominalInterpolate" << std::endl;
            exit(0);
          }
          den = w/den;
          c[i]=(ho*den);
          d[i]=(hp*den);
        }

        fDelta = ((2.0*closestIndex < (dim - mi)) ? c[closestIndex+1] :
                 d[closestIndex--]);
        yOut += fDelta;
      }

      return yOut;
    }

  static G4double MThPolynominalInterpolate(const std::vector<double>& xa, const std::vector<double>& ya,
					    int n,
					    G4double x,
					    int mi)
    {
    //******************************************************************
    // Interpolate with mTh pylinominal function. The data arrays
    // xa[0.....n-1] and ya[0....n-1] must be larger than m,
    // polynominal of degree.
    //******************************************************************

      int index = SearchIndex(xa, x, n);
      int kMax = ((index-(mi-1)/2 > 1) ? index-(mi-1)/2 : 1);
      int kMin = ((kMax < n+1-mi) ? kMax : n+1-mi);

      std::vector<double> xpart;
      std::vector<double> ypart;

      xpart.resize(mi);
      ypart.resize(mi);

      for(int i=0; i<mi; i++){
        xpart[i]=xa[kMin-1+i];
        ypart[i]=ya[kMin-1+i];
      }


      double yOut = PolynominalInterpolate(xpart,ypart,mi,x);
      //std::cout << "+++++++++" << yOut << std::endl;
      return yOut;

    }



  static void DataReader(const G4String inputfname, std::vector<double> &xaxis, std::vector<double> &yaxis)
    {
    //
    // Read inputfile and store data into two std::vector<double>.
    // Note that these std::vector<double> (xaxis, yaxis) are cleard before
    // read data.
    //
       // clear current array; set 0 for all components
       xaxis.clear();
       yaxis.clear();

       std::ifstream in(inputfname);
       if (!in.good()) {
          std::cerr << "Interpolation:DataReader: file " << inputfname
               << " open failed! abort." << std::endl;
          abort();
       }

       const int nchar = 2048;
       char      buf[nchar];
       std::string     dum1;
       double    dumx, dumy;

       int  ndata = -1;  // initialize

       /*while (!in.eof()) {

         if (ndata != -1) {  // ndata is set already
            for (int i=0; i<ndata; i++) {
                //in >> dumx >> dumy;*/
               /* xaxis.push_back(dumx);
                yaxis.push_back(dumy);

                xaxis[i] = dumx;
                yaxis[i] = dumy;
             }
             // now exit while loop
             //break;
          }

          // ndata is not set yet.
          // read ndata and skip headers

          in.getline(buf,nchar);

          dum1 = buf;
          std::string dum2;

          if (dum1.find("ndata") != std::string::npos) {
             std::stringstream ss(buf);
             ss >> dum2 >> ndata;

             if (ndata < 0) {
                std::cerr << "Interpolation:DataReader: "
                     << "invalid ndata! ndata = " << ndata << std::endl;
                abort();
             }

             // set size of array and allocate memory
             xaxis.resize(ndata);
             yaxis.resize(ndata);
             continue;

          } else {
             // Header or LF. do nothing.
          }

       } */

       while(in >> dumx >> dumy)
       {
        xaxis.push_back(dumx);
        yaxis.push_back(dumy);
       }

       in.close();
       in.clear();
    }


  static double GetErrorInPolynominalInterpolate( );

     // error of interpolation

};

#endif
