#include "Interpolation.hh"

double Interpolation::fDelta = 1.0;

//__________________________________________________________________
double Interpolation::GetErrorInPolynominalInterpolate( ) {
  return fDelta;
}

