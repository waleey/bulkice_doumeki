#include "OMSimRadioactivityData.hh"

extern G4bool gVerbose;
/**
*Time window will be set in Particle Setup Class
*Out and In rad will be set in Detector Construction.
**/

//G4double OMSimRadioactivityData::ftimeWindow = 10;
//G4double OMSimRadioactivityData::fglassInRad = 10;
//G4double OMSimRadioactivityData::fglassOutRad = 10;
OMSimRadioactivityData::OMSimRadioactivityData() : factivity(0), fnumDecay(0), finitialTime(0)
{
}
OMSimRadioactivityData::~OMSimRadioactivityData()
{
}

G4ThreeVector OMSimRadioactivityData::SetupPosition()
{
    GeneratePosition();
    return fdecayPosition;
}

G4ThreeVector OMSimRadioactivityData::SetupOrientation()
{
    GenerateOrientation();
    return fdecayOrientation;
}

void OMSimRadioactivityData::SetTimeWindow(G4double timeWindow)
{
    OMSimRadioactivityData::ftimeWindow = timeWindow;
}

void OMSimRadioactivityData::SetGlassRad(G4double inRad, G4double outRad)
{
    OMSimRadioactivityData::fglassOutRad = outRad;
    OMSimRadioactivityData::fglassInRad = inRad;
}

G4int OMSimRadioactivityData::GetNumDecay(G4double meanRate)
{
    InversePoisson(meanRate);
    return int(fnumDecay);
}

G4double OMSimRadioactivityData::GetTimeWindow()
{
    return ftimeWindow;
}

G4double OMSimRadioactivityData::GetInitialTime()
{
    GenerateFlatInitialTime();
    //GenerateExpInitialTime();

    return finitialTime;
}

G4double OMSimRadioactivityData::GetInitialTime(G4String type, G4double meanLifeTime)
{
    if (type == "flat") {GenerateFlatInitialTime();}
    else if (type == "exp") {GenerateExpInitialTime(meanLifeTime);}
    return finitialTime;
}

void OMSimRadioactivityData::GenerateFlatInitialTime()
{
    finitialTime = RandomGen(0.0, ftimeWindow);
}

G4double OMSimRadioactivityData::GetInitialTimeBounds(G4double timeLow, G4double timeHigh)
{
    finitialTime = RandomGen(timeLow, timeHigh);
    return finitialTime;
}

void OMSimRadioactivityData::GenerateExpInitialTime(G4double meanLifeTime)
{
    // inverse CDF sampling for normalized PDF between [0, ftimeWindow] and given meanLifeTime
    // for large fTimeWindow > 10E16 s floating point precision issues emerge
    // switch to flat distribution for those cases
    G4double randomCDF = RandomGen(0.0, 1.0);
    if (gVerbose){
        std::cout << "+++ (RAD): randomCDF = " << randomCDF << std::endl;
    }
    finitialTime = -meanLifeTime * log(1 - randomCDF * (1 - exp(-ftimeWindow/meanLifeTime)));
}

void OMSimRadioactivityData::GeneratePosition()
{
    G4double xIn;
    G4double yIn;
    G4double zIn;

    /**
    *For now, radioactivity is only available to MDOM
    *Other OMs will be added soon once the geometry is well known.
    **/
    G4int omModel = OMSimRadioactivityData::fomModel;
    try
    {
        if( omModel== 1)
        {
            fglassWeight = 13.0; //kilograms
            G4double glassOutRad = OMSimRadioactivityData::fglassOutRad;
            G4double glassInRad = OMSimRadioactivityData::fglassInRad;
            fglassRadius = RandomGen(glassInRad, glassOutRad);
            fglassTheta = RandomGen(0, CLHEP::pi);
            fglassPhi = RandomGen(0, CLHEP::pi * 2);

            xIn = fglassRadius * sin(fglassTheta) * cos(fglassPhi);
            yIn = fglassRadius * sin(fglassTheta) * sin(fglassPhi);
            zIn = fglassRadius * cos(fglassTheta);


            fdecayPosition = G4ThreeVector(xIn * mm, yIn * mm, zIn * mm);
        }
        else
        {
            throw omModel;
        }
    }
    catch(G4int omModel)
    {
        std::cout << "Invalid OM Model. Radioactive features not available for this model yet. Aborting..." << std::endl;
    }
}

void OMSimRadioactivityData::GenerateOrientation()
{
    G4double orientationRandom[3] = {RandomGen(-1,1),RandomGen(-1,1),RandomGen(-1,1)};
    G4double mode = std::sqrt(std::pow(orientationRandom[0] , 2) + std::pow(orientationRandom[1] , 2) + std::pow(orientationRandom[2] , 2));
    G4double orientDirection[3] = {orientationRandom[0]/mode, orientationRandom[1]/mode, orientationRandom[2]/mode};

    fdecayOrientation = G4ThreeVector(orientDirection[0], orientDirection[1], orientDirection[2]);
}

void OMSimRadioactivityData::InversePoisson(G4double meanRate)
{
    G4double minLim = 0.0;
    G4double maxLim = 1.0;

    G4double uniform_random = RandomGen(minLim, maxLim);
    G4double mu = meanRate;
    fnumDecay = 0;

    boost::math::poisson_distribution<> dist(mu);
    while(boost::math::cdf(dist, fnumDecay) < uniform_random)
    {
        ++fnumDecay;
    }

}

G4double OMSimRadioactivityData::RandomGen(G4double minLim, G4double maxLim)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minLim, maxLim);

    return dis(gen);
}
