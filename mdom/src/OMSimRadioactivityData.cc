#include "OMSimRadioactivityData.hh"
/**
*Time window will be set in Particle Setup Class
*Out and In rad will be set in Detector Construction.
**/

//G4double OMSimRadioactivityData::ftimeWindow = 10;
//G4double OMSimRadioactivityData::fglassInRad = 10;
//G4double OMSimRadioactivityData::fglassOutRad = 10;
OMSimRadioactivityData::OMSimRadioactivityData(G4int omModel) : fomModel(omModel), factivity(0), fnumDecay(0), finitialTime(0)
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

G4int OMSimRadioactivityData::GetNumDecay()
{
    InversePoisson();
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

void OMSimRadioactivityData::GenerateFlatInitialTime()
{
    finitialTime = RandomGen(0.0, ftimeWindow);
}

void OMSimRadioactivityData::GenerateExpInitialTime()
{

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

    try
    {
        if(fomModel == 1)
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
            throw fomModel;
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

void OMSimRadioactivityData::InversePoisson()
{
    G4double minLim = 0.0;
    G4double maxLim = 1.0;

    G4double uniform_random = RandomGen(minLim, maxLim);
    G4double mu = ftimeWindow * factivity;
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
