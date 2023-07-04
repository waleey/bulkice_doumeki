#include "Poisson.hh"

void Poisson::InversePoisson()
{
    G4double uniform_random = RandomGen();
    mstepSize = 1.0;
    G4double mu = mtimeWindow * mactivity;
    mnumDecay = 0;

    boost::math::poisson_distribution<> dist(mu);
    while(boost::math::cdf(dist, mnumDecay) < uniform_random)
    {
        ++mnumDecay;
    }

}

G4double Poisson::RandomGen()
{
    mminLim = 0.0;
    mmaxLim = 1.0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(mminLim, mmaxLim);

    return dis(gen);
}

G4int Poisson::GetNumDecay()
{
    InversePoisson();

    return int(mnumDecay);
}
