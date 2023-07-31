#include "OMSimPMTResponse.hh"

OMSimPMTResponse::OMSimPMTResponse(std::string& qe_filename) : fqe_filename(qe_filename)
{
    std::cout << "*************PMT Response class created! QE is considered****************" << std::endl;
    SetupQE();
}
OMSimPMTResponse::~OMSimPMTResponse()
{}

bool OMSimPMTResponse::PhotonRegistered(double wavelen)
{
    double rand_num = RandomGen();
    return ( rand_num < qe_map[wavelen] ) ? true : false;
}
void OMSimPMTResponse::SetupQE()
{
    std::ifstream in_file(fqe_filename);

    if(!in_file)
    {
        std::cerr << "XXXX Failed to open " << fqe_filename << " file XXXX" << std::endl;
    }
    else
    {
        double temp_wavelen, temp_qe;
        while(in_file >> temp_wavelen >> temp_qe)
        {
                qe_map[temp_wavelen] = temp_qe;
        }

        std::cout << "*******QE Data Read Successfully**********" << std::endl;
    }
}

double OMSimPMTResponse::RandomGen()
{
	std::random_device rd;  // Obtain a random seed from the hardware
    std::mt19937 gen(rd()); // Seed the generator

    std::uniform_real_distribution<double> dis(0.0, 1.0); // Define the distribution

    double random_num = dis(gen); // Generate a random number

    //std::cout << "Random number between 0 and 1: " << randomNum << std::endl;

    return random_num;
}
void OMSimPMTResponse::PrintQE() const
{
    std::cout << "PMT QE TABLE" << std::endl;
    for(const auto pair : qe_map)
    {
        std::cout << "Wlen: " << pair.first << " " << "QE" << pair.second << std::endl;
    }
}
