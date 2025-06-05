#define _USE_MATH_DEFINES
#include <cmath>

#include "OMSimNeutrinoAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include "CLHEP/Random/RandGamma.h"
#include "Randomize.hh"

extern G4double gworldsize;
extern G4bool   gVerbose;
extern G4int    gRunID;

extern std::fstream gPositronData;


extern G4double gNeutrinoMeanEnergy;
extern G4double gNeutrinoEnergyPinch;
extern G4String gNeutrinoEnergyType;
extern G4double gPositronDensity;
extern G4String gPositronNumber;
extern G4int    gPositronEnergyOrder;
extern G4int    gPositronZenithOrder;
extern G4String gPositronZenithType;


OMSimNeutrinoAction::OMSimNeutrinoAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0),
        fPerc(0)
{
    fRadData = new OMSimRadioactivityData();
}

OMSimNeutrinoAction::~OMSimNeutrinoAction()
{

}

void OMSimNeutrinoAction::GeneratePrimaries(G4Event* anEvent)
{
    // pull neutrino energy from distribution
    G4double neutrinoEnergy = 0.0 * MeV;
    if (gNeutrinoEnergyType=="mono") // mono-energetic
    {
        neutrinoEnergy = gNeutrinoMeanEnergy;
    }
    else // gamma distribution
    {
        do {
        neutrinoEnergy = SampleEnergy(gNeutrinoMeanEnergy, gNeutrinoEnergyPinch);
        } while (neutrinoEnergy < IBDThresholdEnergy);
    }

    // zeroth order positron energy E_pos = E_nu - (m_neu - m_pro)
    G4double positronEnergy0 = neutrinoEnergy - (neutronMass - protonMass);

    // zeroth order positron momentum
    G4double positronMomentum0 = std::sqrt(positronEnergy0*positronEnergy0 - electronMass*electronMass); 
    
    // zeroth order velocity in units of c
    G4double positronVe0 = positronMomentum0 / positronEnergy0; 

    // pull cos(zenith) from distribution
    G4double positronCosZenith = 0.0;
    if (gPositronZenithType=="uniform") // uniform
    {
        positronCosZenith = fRadData -> RandomGen(-1.0, 1.0);
    }
    else // realistic (asymmetrical) distribution
    {
        positronCosZenith = SampleZenith(positronVe0, neutrinoEnergy, "analytic", gPositronZenithOrder);
    }

    // debug
    G4double positronCosZenithUni = fRadData -> RandomGen(-1.0, 1.0);
    G4double positronCosZenith0 = SampleZenith(positronVe0, neutrinoEnergy, "analytic", 0);
    G4double positronCosZenith1 = SampleZenith(positronVe0, neutrinoEnergy, "analytic", 1);


    // first order positron energy
    G4double positronEnergy1 = positronEnergy0*(1-(neutrinoEnergy/protonMass)*(1-positronVe0*positronCosZenith))-(pow((neutronMass-protonMass),2)-pow(electronMass,2))/(2*protonMass);

    // return positron energy based on order
    G4double positronEnergy = 0.0 * MeV;
    if (gPositronEnergyOrder == 0)
    {
        positronEnergy = positronEnergy0;
    }
    else if (gPositronEnergyOrder == 1)
    {
        positronEnergy = positronEnergy1;   
    }
    
    // positron sin zenith
    G4double positronSinZenith = std::sqrt(1.0 - positronCosZenith * positronCosZenith);

    // downgoing neutrinos
    G4ThreeVector neutrinoDirection(0,0,-1);

    // positron azimuth angle 
    G4double positronAzimuth = fRadData -> RandomGen(0.0, 2*M_PI);

    // positron direction in neutrino frame
    G4ThreeVector positronDirectionInNeutrinoFrame (
        positronSinZenith * std::cos(positronAzimuth),
        positronSinZenith * std::sin(positronAzimuth),
        positronCosZenith
    );

    // rotate positron direction into lab frame, custom method
    G4ThreeVector positronDirection = rotateFrameFromeTo(positronDirectionInNeutrinoFrame, G4ThreeVector(0,0,1), neutrinoDirection);

    // random positron position inside simulation volume
    // for box volume
    //G4ThreeVector positronPosition (
    //    fRadData -> RandomGen(-gworldsize, gworldsize) * m,
    //    fRadData -> RandomGen(-gworldsize, gworldsize) * m,
    //    fRadData -> RandomGen(-gworldsize, gworldsize) * m
    //);
    // for spherical volume
    G4double positronPositionAzimuth = 2*M_PI * G4UniformRand();
    G4double positronPositionZenith = std::acos(1 - 2 * G4UniformRand());
    G4double positronPositionRadius = gworldsize * std::cbrt(G4UniformRand());

    G4ThreeVector positronPosition (
        positronPositionRadius * std::sin(positronPositionZenith) * std::cos(positronPositionAzimuth) * m,
        positronPositionRadius * std::sin(positronPositionZenith) * std::sin(positronPositionAzimuth) * m,
        positronPositionRadius * std::cos(positronPositionZenith) * m
    );

    if (gVerbose)
    {
        std::cout << "Neutrino energy: " << neutrinoEnergy / MeV << " MeV" << "\n";
        std::cout << "Positron energy: " << positronEnergy / MeV << " MeV" << "\n";
        std::cout << "Positron zenith angle: " << std::acos(positronCosZenith)/M_PI * 180 << " deg \n";
        std::cout << "Positron azimuth angle: " << positronAzimuth/M_PI * 180 << " deg \n";
        std::cout << "Neutrino direction: " << neutrinoDirection << "\n";
        std::cout << "Positron direction (neutrino frame): " << positronDirectionInNeutrinoFrame << "\n";
        std::cout << "Positron direction (lab frame): " << positronDirection << "\n";
        std::cout << "Positron position (lab frame): " << positronPosition / m << "\n";
    }

    // save positron energy, cos(zenith) and azimuth to file
    if (gPositronData.is_open())
    {
        gPositronData << neutrinoEnergy / MeV << "," << positronEnergy0 / MeV << "," << positronEnergy1 / MeV << "," 
                      << positronCosZenithUni << "," << positronCosZenith0 << "," << positronCosZenith1 << "," 
                      << positronAzimuth << "," << positronPosition[0] / m << "," << positronPosition[1] / m << "," 
                      << positronPosition[2] / m << "\n";
        gPositronData.flush();
    }
    
    // sample time from uniform distribution
    G4double timeRandom = fRadData -> RandomGen(0, 1);

    // set variables for particle gun
    fParticleGun -> SetParticlePosition(positronPosition);
    fParticleGun -> SetParticleMomentumDirection(positronDirection);
    fParticleGun -> SetParticleEnergy(positronEnergy);
    fParticleGun -> SetParticleTime(timeRandom * s);
    fParticleGun -> SetParticleDefinition(G4Positron::PositronDefinition());
    fParticleGun -> GeneratePrimaryVertex(anEvent);

    // prints progress for every percentage point
    if (fIdx > (fParticleNum * fPerc/100)){
        std::cerr << fPerc << " percent done" << std::endl;
        fPerc++;
    }

    fIdx++;
    if(fIdx == fParticleNum)
    {
        std::cout << "Total " << fIdx << " positrons are generated!" << std::endl;
        fParticleExist = false;
    }
}

G4double OMSimNeutrinoAction::SampleEnergy(G4double meanEnergy, G4double alpha)
{
    G4double k = alpha + 1;
    G4double theta = meanEnergy / k;
    G4double lambda = 1/theta;
    // CLHEP::RandGamma uses k and lambda = 1/theta
    return CLHEP::RandGamma::shoot(k, lambda);; // uses global RNG
}

G4double OMSimNeutrinoAction::SampleZenith(G4double positronVe0, G4double neutrinoEnergy, G4String mode, G4int order)
{  
    // mean of the cos(zenith) distribution
    G4double mu0 = 1.0/3.0 * positronVe0 * a0; //zeroth order
    G4double mu1 = mu0 + 2.4 * neutrinoEnergy/protonMass; // first order
    G4double mu = mu0;
    if (order == 1){mu = mu1;}

    // pdf(cos(zenith)) = 1/2 ( 1 + 3 µ cos(zenith))
    // cdf(x) = int_-1^x(pdf) = 1/2 ((x+1) + 3/2 µ (x*x -1)) = u

    if (mode=="analytic")
    {
        // numeric, inverse CDF
        // rearrange cdf(x) in terms of x

        G4double A = 0.75 * mu;
        G4double B = 0.5;
        G4double C = 0.5 - 0.75 * mu;

        G4double U = fRadData -> RandomGen(0.0, 1.0); // random CDF value between 0 and 1

        // solve quadratic equation
        G4double cos1 = - (B/(2*A)) + std::sqrt(pow(B/(2*A), 2) - (C-U)/A);
        G4double cos2 = - (B/(2*A)) - std::sqrt(pow(B/(2*A), 2) - (C-U)/A);

        // one root is unphysical
        if (std::abs(cos1) > 1){ return cos2; }
        else { return cos1; }
    }
    else
    {
        G4double pdf_max = 0.5 * (1 + 3 * std::abs(mu)); // maximum of the pdf
        G4double x;
        G4double y;

        while (true) 
        {
            x = fRadData -> RandomGen(-1.0, 1.0);
            y = fRadData -> RandomGen(0, pdf_max);
            double f = 0.5 * (1 + 3 * mu * x);
            if (y <= f) return x;
        }
    }
}

G4ThreeVector OMSimNeutrinoAction::rotateFrameFromeTo(const G4ThreeVector& original, const G4ThreeVector& from, const G4ThreeVector& to)
{
    // normalize from and to vector
    G4ThreeVector v = from.unit();
    G4ThreeVector w = to.unit();

    // get orthogonal rotation axis
    G4ThreeVector axis = v.cross(w);
    G4double sinAngle = axis.mag();
    G4double cosAngle = v.dot(w);

    // degenerate case: to and from are parallel/anti-parallel
    if (sinAngle < 1e-12)
    {
        if (cosAngle > 0.999999) 
        {
            return original; // parallel, no flip
        }
        else
        {
            G4ThreeVector flipped = G4ThreeVector(-original.x(), original.y(), -original.z()); // rotation around y-axis
            return flipped;
        }
    }
    // general case: perform rotation
    axis = axis.unit();
    G4double angle = std::atan2(sinAngle, cosAngle);
    G4ThreeVector rotated = original;
    rotated.rotate(angle, axis);
    return rotated;
}

void OMSimNeutrinoAction::LoadData()
{
    G4double meanRate = gPositronDensity * pow((gworldsize * 2),3);

    if (gPositronNumber == "poisson") { fParticleNum = fRadData -> GetNumDecay(meanRate); }
    else {fParticleNum = meanRate; }

    std::cout << "---------------------\n" 
              << "Positron density [pos/m3]: " << gPositronDensity << "\n"
              << "Simulation Volume [m]: " << pow((gworldsize * 2),3) << "\n"
              << "Number of expected positrons: " << meanRate << "\n"
              << "Number of simulated positrons: " << fParticleNum << "\n"
              << "---------------------" << std::endl;

}