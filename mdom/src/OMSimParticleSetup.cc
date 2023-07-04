#include "OMSimParticleSetup.hh"

OMSimParticleSetup::OMSimParticleSetup(G4ParticleGun* ParticleGun, G4Event* anEvent, G4int omModel) : fParticleGun(ParticleGun), fEvent(anEvent), fomModel(omModel), fglassweight(0), timeWindow(60)
{
    /**
    *For now, time winodw is set to 1 min.
    *Should be changed later to 10 min as M.Unland's thesis.
    **/
    poisson = new Poisson(1, 1);
}
void OMSimParticleSetup::GeneratePositron()
{
    SetupPositron();
    G4String particle_name = "e+";
    G4ParticleDefinition* particle = fParticleTable -> FindParticle(particle_name);
    int pc{0};
    for(G4int i = 0; i < pos_count; i++)
	{
		G4ThreeVector particlePosition(pos_data[X][i] * m, pos_data[Y][i] * m, pos_data[Z][i] * m);
		G4ThreeVector particleOrientation(pos_data.at(AX)[i], pos_data.at(AY)[i], pos_data.at(AZ)[i]);

		G4double particleEnergy = pos_data[ENERGY][i] * CLHEP::MeV;
		G4double particleInTime = pos_data[TIME][i] *ms;
		fParticleGun -> SetParticlePosition(particlePosition);
		fParticleGun -> SetParticleMomentumDirection(particleOrientation);
		fParticleGun -> SetParticleEnergy(particleEnergy);
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleTime(particleInTime);
		fParticleGun -> GeneratePrimaryVertex(fEvent);
		pc++;
	}
	std::cout << "Total PC: " << pc << std::endl;
	//clear the data vector and get ready for next run
	for(auto& vals: pos_data)
	{
        vals.clear();
	}

}

void OMSimParticleSetup::GenerateNeutron()
{
    SetupNeutron();
    G4String particle_name = "neutron";
    G4ParticleDefinition* particle = fParticleTable -> FindParticle(particle_name);
    for(G4int i = 0; i < neu_count; i++)
	{
		G4ThreeVector particlePosition(neu_data[X][i] * m, neu_data[Y][i] * m, neu_data[Z][i] * m);
		G4ThreeVector particleOrientation(neu_data.at(AX)[i], neu_data.at(AY)[i], neu_data.at(AZ)[i]);
		G4double particleEnergy = neu_data[ENERGY][i] * CLHEP::MeV;
		G4double particleInTime = neu_data[TIME][i] *ms;

		fParticleGun -> SetParticlePosition(particlePosition);
		fParticleGun -> SetParticleMomentumDirection(particleOrientation);
		fParticleGun -> SetParticleEnergy(particleEnergy);
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleTime(particleInTime);
		fParticleGun -> GeneratePrimaryVertex(fEvent);

	}
	//clear the vector and get ready for next run
	for(auto& vals: neu_data)
	{
        vals.clear();
	}

}

void OMSimParticleSetup::GenerateElectron()
{
    SetupElectron();

    G4String particle_name = "e-";
    G4ParticleDefinition* particle = fParticleTable -> FindParticle(particle_name);
    for(G4int i = 0; i < e_count; i++)
	{
		G4ThreeVector particlePosition(e_data[X][i] * m, e_data[Y][i] * m, e_data[Z][i] * m);
		G4ThreeVector particleOrientation(e_data.at(AX)[i], e_data.at(AY)[i], e_data.at(AZ)[i]);

		G4double particleEnergy = e_data[ENERGY][i] * CLHEP::MeV;
		G4double particleInTime = e_data[TIME][i] *ms;
		fParticleGun -> SetParticlePosition(particlePosition);
		fParticleGun -> SetParticleMomentumDirection(particleOrientation);
		fParticleGun -> SetParticleEnergy(particleEnergy);
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleTime(particleInTime);
		fParticleGun -> GeneratePrimaryVertex(fEvent);

	}
	//clear the vector and get ready for next run
	for(auto& vals: e_data)
	{
        vals.clear();
	}

}
void OMSimParticleSetup::GenerateK40()
{
    G4int Z = 19;
    G4int A = 40;
    SetupRadioactiveParticles();
    k40activity = 61; /** value from C. Lozano's 2016 paper**/
    poisson -> SetTimeWindow(timeWindow);
    poisson -> SetActivity(k40activity);
    radioactiveParticleNum = poisson -> GetNumDecay();
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);

    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        SetupRadioactiveParticles();
        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}
void OMSimParticleSetup::GenerateTh238()
{
    G4int Z = 90;
    G4int A = 238;
    SetupRadioactiveParticles();
    th238activity = 1.28; /** value from C. Lozano's 2016 paper**/
    poisson -> SetTimeWindow(timeWindow);
    poisson -> SetActivity(th238activity);
    radioactiveParticleNum = poisson -> GetNumDecay();
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);
    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        SetupRadioactiveParticles();
        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}
void OMSimParticleSetup::GenerateU238()
{
    G4int Z = 92;
    G4int A = 238;
    SetupRadioactiveParticles();
    u238activity = 4.61; /** value from C. Lozano's 2016 paper **/
    poisson -> SetTimeWindow(timeWindow);
    poisson -> SetActivity(u238activity);
    radioactiveParticleNum = poisson -> GetNumDecay();
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);
    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        SetupRadioactiveParticles();
        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}
void OMSimParticleSetup::GenerateU235()
{
    G4int Z = 92;
    G4int A = 235;
    SetupRadioactiveParticles();
    u235activity = 0.59; /** value from C. Lozano's 2016 paper **/
    poisson -> SetTimeWindow(timeWindow);
    poisson -> SetActivity(u235activity);
    radioactiveParticleNum = poisson -> GetNumDecay();
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);
    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        SetupRadioactiveParticles();
        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}

void OMSimParticleSetup::SetupPositron()
{
    DataReader(pos_filePath, pos_data);
    pos_count = pos_data.at(ENERGY).size();
    std::cout << "Positron information are set up for " << pos_count << " positrons!" << std::endl;

}

void OMSimParticleSetup::SetupNeutron()
{
    DataReader(neu_filePath, neu_data);
    neu_count = neu_data.at(ENERGY).size();
    std::cout << "Neutron Information are set up for " << neu_count << " neutrons!" << std::endl;
}

void OMSimParticleSetup::SetupElectron()
{
    DataReader(e_filePath, e_data);
    e_count = e_data.at(ENERGY).size();
    std::cout << "Electron Information are set up for " << e_count << " electrons!" <<std::endl;
}
void OMSimParticleSetup::SetupRadioactiveParticles()
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
            fglassweight = 13.0; //kilograms
            G4double glassOutRad = OMSimParticleSetup::glassOutRad;
            G4double glassInRad = OMSimParticleSetup::glassInRad;
            G4double radius = RandomGen(glassInRad, glassOutRad);
            G4double theta = RandomGen(0, CLHEP::pi);
            G4double phi = RandomGen(0, CLHEP::pi * 2);

            xIn = radius * sin(theta) * cos(phi);
            yIn = radius * sin(theta) * sin(phi);
            zIn = radius * cos(theta);

            G4double orientationRandom[3] = {RandomGen(-1,1),RandomGen(-1,1),RandomGen(-1,1)};
            G4double mode = std::sqrt(std::pow(orientationRandom[0],2)+std::pow(orientationRandom[1],2)+std::pow(orientationRandom[2],2));
            G4double orientDirection[3] = {orientationRandom[0]/mode, orientationRandom[1]/mode, orientationRandom[2]/mode};

            radioactiveParticleEnergy = 0*keV;
            radioactiveParticleCharge = 0.*eplus;
            radioactiveParticlePosition = G4ThreeVector(xIn, yIn, zIn);
            radioactiveParticleOrientation = G4ThreeVector(orientDirection[0], orientDirection[1], orientDirection[2]);
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
void OMSimParticleSetup::DataReader(std::string filePath, std::vector<std::vector<G4double>>& data)
{
    using namespace std;
    G4double temp;
    string fileName;

    for(int i = 0; i < data.size(); i++)
    {
        temp = 0;
        fileName = filePath + dtypes.at(i) + ".data";

        ifstream file(fileName);

	if(!file.is_open())
	{
		std::cout << "Failed to open" + fileName << std::endl;

	}

	while(file >> temp)
	{
        data.at(i).push_back(temp);
	}

	file.close();
    }

}
double OMSimParticleSetup::RandomGen(double minLim, double maxLim)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minLim, maxLim);

    return dis(gen);
}
void OMSimParticleSetup::SetGlassRad(G4double outrad, G4double inrad)
{
    OMSimParticleSetup::glassInRad = inrad;
    OMSimParticleSetup::glassOutRad = outrad;
}
