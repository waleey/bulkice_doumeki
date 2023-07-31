#ifndef OMSIMPMTRESPONSE_HH_INCLUDED
#define OMSIMPMTRESPONSE_HH_INCLUDED

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <random>

class OMSimPMTResponse
{
public:
    OMSimPMTResponse(std::string& qe_filename);
    ~OMSimPMTResponse();

    double GetQE() const;
    bool PhotonRegistered(double wavelen);
    void PrintQE() const;

private:
    double fwavelen;
    double fqe;
    std::string fqe_filename;
    std::map<double, double> qe_map;

    void SetupQE();
    double RandomGen();

};
#endif

