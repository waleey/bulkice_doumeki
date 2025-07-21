#define _USE_MATH_DEFINES
#include <cmath>

#include "OMSimNeutrinoAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include "CLHEP/Random/RandGamma.h"
#include "Randomize.hh"
#include "spline.h"

extern G4double gworldsize;
extern G4bool   gVerbose;
extern G4int    gRunID;

//extern std::fstream gPositronData;

extern G4int    gCrossSectionOrder;
extern G4int    gPositronEnergyOrder;
extern G4int    gPositronZenithOrder;

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
    G4double neutrinoEnergy = SampleNeutrinoEnergy();

    // zeroth order positron energy and momentum
    G4double positronEnergyZeroth = IBDPositronEnergyZeroth(neutrinoEnergy);
    G4double positronMomentumZeroth = IBDPositronMomentumZeroth(positronEnergyZeroth);

    // pull cos(zenith) from distribution
    G4double positronCosZenith = SamplePositronZenith(neutrinoEnergy, positronEnergyZeroth, "analytic", gPositronZenithOrder);

    // debug
    G4double positronCosZenithUni = fRadData -> RandomGen(-1.0, 1.0);
    G4double positronCosZenithZeroth = SamplePositronZenith(neutrinoEnergy, positronEnergyZeroth, "analytic", 0);
    G4double positronCosZenithFirst = SamplePositronZenith(neutrinoEnergy, positronEnergyZeroth, "analytic", 1);

    // first order positron energy
    G4double positronEnergyFirst = IBDPositronEnergyFirst(neutrinoEnergy, positronEnergyZeroth, positronCosZenith);

    // return positron energy based on order
    G4double positronEnergy = 0.0 * MeV;
    if (gPositronEnergyOrder == 0)
    {
        positronEnergy = positronEnergyZeroth;
    }
    else if (gPositronEnergyOrder == 1)
    {
        positronEnergy = positronEnergyFirst;   
    }
    
    // positron sin zenith
    G4double positronSinZenith = std::sqrt(1.0 - positronCosZenith * positronCosZenith);

    // downgoing neutrinos
    G4ThreeVector neutrinoDirection(0,0,-1);

    // positron azimuth angle 
    G4double positronAzimuth = CLHEP::RandFlat::shoot(0.0,1.0);

    // positron direction in neutrino frame
    G4ThreeVector positronDirectionInNeutrinoFrame (
        positronSinZenith * std::cos(positronAzimuth),
        positronSinZenith * std::sin(positronAzimuth),
        positronCosZenith
    );

    // rotate positron direction into lab frame, custom method
    G4ThreeVector positronDirection = rotateFrameFromeTo(positronDirectionInNeutrinoFrame, G4ThreeVector(0,0,1), neutrinoDirection);

    // random positron position inside simulation volume
    // for spherical volume
    G4double positronPositionAzimuth = 2*M_PI * G4UniformRand();
    G4double positronPositionZenith = std::acos(1 - 2 * G4UniformRand());
    G4double positronPositionRadius = fworldsize * std::cbrt(G4UniformRand());

    G4ThreeVector positronPosition (
        positronPositionRadius * std::sin(positronPositionZenith) * std::cos(positronPositionAzimuth),
        positronPositionRadius * std::sin(positronPositionZenith) * std::sin(positronPositionAzimuth),
        positronPositionRadius * std::cos(positronPositionZenith)
    );

    if (gVerbose)
    {
        std::cout << "Neutrino energy [MeV]: " << neutrinoEnergy / MeV << " MeV" << "\n";
        std::cout << "Positron energy [MeV]: " << positronEnergy / MeV << " MeV" << "\n";
        std::cout << "Positron zenith angle [rad]: " << std::acos(positronCosZenith)/M_PI * 180 << " deg \n";
        std::cout << "Positron azimuth angle [rad]: " << positronAzimuth/M_PI * 180 << " deg \n";
        std::cout << "Neutrino direction: " << neutrinoDirection << "\n";
        std::cout << "Positron direction (neutrino frame): " << positronDirectionInNeutrinoFrame << "\n";
        std::cout << "Positron direction (lab frame): " << positronDirection << "\n";
        std::cout << "Positron position (lab frame) [m]: " << positronPosition / m << "\n";
    }

    // save positron energy, cos(zenith) and azimuth to file
    //if (gPositronData.is_open())
    //{
    //    gPositronData << neutrinoEnergy / MeV << "," << positronEnergyZeroth / MeV << "," << positronEnergyFirst / MeV << "," 
    //                  << positronCosZenithUni << "," << positronCosZenithZeroth << "," << positronCosZenithFirst << "," 
    //                  << positronAzimuth << "," << positronPosition[0] / m << "," << positronPosition[1] / m << "," 
    //                  << positronPosition[2] / m << "\n";
    //    gPositronData.flush();
    //}
    
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

G4double OMSimNeutrinoAction::SampleNeutrinoEnergy()
{
    G4double u = CLHEP::RandFlat::shoot(0.0,1.0);
    return fInverseCDF(u); // interpolated energy
}

G4double OMSimNeutrinoAction::SamplePositronZenith(G4double neutrinoEnergy, G4double positronEnergyZeroth, G4String mode, G4int order)
{  
    G4double positronVeZeroth = IBDPositronVeZeroth(positronEnergyZeroth);
    // mean of the cos(zenith) distribution
    G4double mu0 = 1.0/3.0 * positronVeZeroth * a0; //zeroth order
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

void OMSimNeutrinoAction::LoadData()
{
    
    std::cout << "---------------------\n"
              << "Neutrino Flux [J/s/cm2] ([MeV/ns/mm2]): " << fNeutrinoFlux / (joule/(s * cm * cm)) 
              << " (" << fNeutrinoFlux / (MeV/(ns * mm * mm)) << ")" <<  "\n"
              << "Neutrino Mean Energy [MeV]: " << fNeutrinoMeanEnergy << "\n"
              << "Neutrino Energy Pinch: " << fNeutrinoEnergyPinch << "\n"
              << "---------------------" << std::endl;
    
    G4double volume = SimulationVolume();
    auto [meanNumberInteraction, numberInteraction] = NumberNeutrinoInteractions();
    InverseCDFNeutrinoSpectrumCrossSection();

    std::cout << "---------------------\n" 
              << "Positron density [pos/m3]: " << numberInteraction/volume * (m*m*m) << "\n"
              << "Simulation Volume [m]: " << volume / (m*m*m) << "\n"
              << "Number of expected positrons: " << meanNumberInteraction << "\n"
              << "Number of simulated positrons: " << numberInteraction << "\n"
              << "---------------------" << std::endl;

    fParticleNum = numberInteraction;
}

void OMSimNeutrinoAction::InverseCDFNeutrinoSpectrumCrossSection()
{
    auto result = NormalizedPDFNeutrinoSpectrumCrossSection(fNeutrinoMeanEnergy, fNeutrinoEnergyPinch, gCrossSectionOrder);

    std::vector<double> energy = result.first;
    std::vector<double> pdf = result.second;

    std::vector<double> cdf = CDFNeutrinoSpectrumCrossSection(pdf);

    fInverseCDF.set_points(cdf, energy); // construct cubic spline E(u)
}

std::vector<double> OMSimNeutrinoAction::CDFNeutrinoSpectrumCrossSection(const std::vector<double>& pdf)
{
    std::vector<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin());
    return cdf;
}

std::pair<std::vector<double>, std::vector<double>> OMSimNeutrinoAction::NormalizedPDFNeutrinoSpectrumCrossSection(G4double neutrinoMeanEnergy, G4double neutrinoEnergyPinch, G4int order)
{
    G4double Emin = IBDThresholdEnergy - 0.01 * MeV; // avoid zero
    G4double Emax = 5 * neutrinoMeanEnergy;
    G4double dE = 0.1 * MeV; // integration step size

    std::vector<double> pdf;
    std::vector<double> energy;
    G4double norm = 0;

    for (G4double E = Emin; E < Emax; E += dE)
    {
        G4double spectrum = NeutrinoBlackBodySpectrum(E, neutrinoMeanEnergy, neutrinoEnergyPinch); // [1/MeV]
        G4double xsec = IBDCrossSection(E, order); // [cm*cm]

        pdf.push_back(spectrum * xsec * dE); // unnormalized pdf
        energy.push_back(E);
        norm += spectrum * xsec * dE;
    }
    for (auto& val : pdf) 
        val /= norm; // dimensionless, normalized to unit area

    return {energy, pdf}; 
}

std::pair<G4double, G4double> OMSimNeutrinoAction::NumberNeutrinoInteractions()
{
    G4double integral = IntegralNeutrinoSpectrumCrossSection(fNeutrinoMeanEnergy, fNeutrinoEnergyPinch, gCrossSectionOrder);
    G4double targetNumber = TargetNumber();

    G4double meanNumberInteraction = fSimulationTime * targetNumber * fNeutrinoFlux / fNeutrinoMeanEnergy * integral; // mean number of neutrino interactions

    G4int numberInteraction = CLHEP::RandPoisson::shoot(meanNumberInteraction); // "real" number of neutrino interactions needs to be pulled from Poisson distribution
    return {meanNumberInteraction, numberInteraction};
}

G4double OMSimNeutrinoAction::IntegralNeutrinoSpectrumCrossSection(G4double neutrinoMeanEnergy, G4double neutrinoEnergyPinch, G4int order)
{
    G4double Emin = 0.01 * MeV; // avoid zero
    G4double Emax = 5 * neutrinoMeanEnergy;
    G4double dE = 0.1 * MeV; // integration step size

    G4double integral = 0.0;

    for (G4double E = Emin; E < Emax; E += dE)
    {
        G4double spectrum = NeutrinoBlackBodySpectrum(E, neutrinoMeanEnergy, neutrinoEnergyPinch); // [1/MeV]
        G4double xsec = IBDCrossSection(E, order); // [cm*cm]
        integral += spectrum * xsec * dE; // [cm*cm]
    }
    return integral;
}

G4double OMSimNeutrinoAction::NeutrinoBlackBodySpectrum(G4double neutrinoEnergy, G4double neutrinoMeanEnergy, G4double neutrinoEnergyPinch)
{
    if (neutrinoEnergy < 0) return 0.0; // guarantees numerical stability

    G4double k = 1 + neutrinoEnergyPinch;
    G4double theta = neutrinoMeanEnergy/k;

    G4double pdf = 1/(std::tgamma(k) * std::pow(theta, k)) * std::pow(neutrinoEnergy, k-1) * std::exp(-neutrinoEnergy/theta);
    return pdf;
}

G4double OMSimNeutrinoAction::IBDCrossSection(G4double neutrinoEnergy, G4int order)
{
    if (order == 0) return IBDCrossSectionZeroth(neutrinoEnergy);
    else return IBDCrossSectionFirst(neutrinoEnergy);
}

G4double OMSimNeutrinoAction::IBDCrossSectionZeroth(G4double neutrinoEnergy)
{
    // From Vogel and Beacom [1999]
    if (neutrinoEnergy < IBDThresholdEnergy) return 0.0;

    G4double positronEnergyZeroth = IBDPositronEnergyZeroth(neutrinoEnergy);
    if (positronEnergyZeroth < electronMass) return 0.0; // guarantees numerical stability

    G4double positronMomentumZeroth = IBDPositronMomentumZeroth(positronEnergyZeroth);

    return 0.0952 * 1E-42 * (positronEnergyZeroth * positronMomentumZeroth)/(MeV * MeV) * (cm*cm);
}

G4double OMSimNeutrinoAction::IBDCrossSectionFirst(G4double neutrinoEnergy)
{
    // From Strumia and Vissani [2003]
    if (neutrinoEnergy < IBDThresholdEnergy) return 0.0;

    G4double positronEnergyZeroth = IBDPositronEnergyZeroth(neutrinoEnergy);
    if (positronEnergyZeroth < electronMass) return 0.0; // guarantees numerical stability

    G4double positronMomentumZeroth = IBDPositronMomentumZeroth(positronEnergyZeroth);

    G4double a = -0.07056;
    G4double b = 0.02018;
    G4double c = -0.001953;

    G4double lnE = std::log(neutrinoEnergy / MeV);
    G4double powerE = std::pow(neutrinoEnergy / MeV, a + b * lnE + c * std::pow(lnE, 3));

    return 1E-43 * (positronEnergyZeroth * positronMomentumZeroth * powerE)/(MeV * MeV) * (cm*cm);
}

G4double OMSimNeutrinoAction::IBDPositronEnergyZeroth(G4double neutrinoEnergy)
{
    return neutrinoEnergy - (neutronMass - protonMass);
}

G4double OMSimNeutrinoAction::IBDPositronMomentumZeroth(G4double positronEnergyZeroth)
{
    return std::sqrt(positronEnergyZeroth*positronEnergyZeroth - electronMass*electronMass);
}

G4double OMSimNeutrinoAction::IBDPositronVeZeroth(G4double positronEnergyZeroth)
{
    G4double positronMomentumZeroth = IBDPositronMomentumZeroth(positronEnergyZeroth);
    return positronMomentumZeroth / positronEnergyZeroth;
}

G4double OMSimNeutrinoAction::IBDPositronEnergyFirst(G4double neutrinoEnergy, G4double positronEnergyZeroth, G4double positronCosZenith)
{
    G4double positronVeZeroth =IBDPositronVeZeroth(positronEnergyZeroth);
    return positronEnergyZeroth * (1-(neutrinoEnergy/protonMass)*(1-positronVeZeroth*positronCosZenith))-(pow((neutronMass-protonMass),2)-pow(electronMass,2))/(2*protonMass);
}

G4double OMSimNeutrinoAction::SimulationVolume()
{
    return (4*M_PI/3 * std::pow(fworldsize, 3));
}

G4double OMSimNeutrinoAction::TargetNumber()
{
    G4double simulationVolume = SimulationVolume();
    G4double targetDensity = 2 * iceDensity/waterMolarMass * avogadroNumber;
    return simulationVolume * targetDensity; // [unitless]
}