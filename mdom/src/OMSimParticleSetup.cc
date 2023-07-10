#include "OMSimParticleSetup.hh"
/**
*Set the time window here
*Do not assign any unit
*It's in seconds
**/
G4double OMSimRadioactivityData::ftimeWindow = 60;

OMSimParticleSetup::OMSimParticleSetup(G4ParticleGun* ParticleGun, G4Event* anEvent, G4int omModel) : fParticleGun(ParticleGun), fEvent(anEvent), fomModel(omModel)
{
    /**
    *For now, time winodw is set to 1 min.
    *Should be changed later to 10 min as M.Unland's thesis.
    **/
    //poisson = new Poisson(1, 1);
    radData = new OMSimRadioactivityData();
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
    k40activity = 61; /** value from C. Lozano's 2016 paper**/
    radData -> SetTimeWindow(OMSimRadioactivityData::ftimeWindow);
    radData -> SetActivity(k40activity);
    radioactiveParticleNum = radData -> GetNumDecay();
    radioactiveParticleEnergy = 0*keV;
    radioactiveParticleCharge = 0.*eplus;
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);

    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        radioactiveParticlePosition = radData -> SetupPosition();
        radioactiveParticleOrientation = radData -> SetupOrientation();
        radioactiveParticleTime = radData -> GetInitialTime() * s;

        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> SetParticleTime(radioactiveParticleTime);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}
void OMSimParticleSetup::GenerateTh238()
{
    G4int Z = 90;
    G4int A = 232;
    th238activity = 1.28; /** value from C. Lozano's 2016 paper**/
    radData -> SetTimeWindow(OMSimRadioactivityData::ftimeWindow);
    radData -> SetActivity(th238activity);
    radioactiveParticleNum = radData -> GetNumDecay();
    radioactiveParticleEnergy = 0*keV;
    radioactiveParticleCharge = 0.*eplus;
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);
    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        radioactiveParticlePosition = radData -> SetupPosition();
        radioactiveParticleOrientation = radData -> SetupOrientation();
        radioactiveParticleTime = radData -> GetInitialTime() * s;

        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> SetParticleTime(radioactiveParticleTime);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}
void OMSimParticleSetup::GenerateU238()
{
    G4int Z = 92;
    G4int A = 238;
    u238activity = 4.61; /** value from C. Lozano's 2016 paper **/
    radData -> SetTimeWindow(OMSimRadioactivityData::ftimeWindow);
    radData -> SetActivity(u238activity);
    radioactiveParticleNum = radData -> GetNumDecay();
    radioactiveParticleEnergy = 0*keV;
    radioactiveParticleCharge = 0.*eplus;
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);
    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        radioactiveParticlePosition = radData -> SetupPosition();
        radioactiveParticleOrientation = radData -> SetupOrientation();
        radioactiveParticleTime = radData -> GetInitialTime() * s;

        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> SetParticleTime(radioactiveParticleTime);
        fParticleGun -> GeneratePrimaryVertex(fEvent);
    }
}
void OMSimParticleSetup::GenerateU235()
{
    G4int Z = 92;
    G4int A = 235;
    u235activity = 0.59; /** value from C. Lozano's 2016 paper **/
    radData -> SetTimeWindow(OMSimRadioactivityData::ftimeWindow);
    radData -> SetActivity(u235activity);
    radioactiveParticleNum = radData -> GetNumDecay();
    radioactiveParticleEnergy = 0*keV;
    radioactiveParticleCharge = 0.*eplus;
    radioactiveParticle = G4IonTable::GetIonTable() -> GetIon(Z, A, radioactiveParticleEnergy);
    radioactiveParticle -> SetPDGLifeTime(0 * ns);
    for(int i = 0; i < radioactiveParticleNum; i++)
    {
        radioactiveParticlePosition = radData -> SetupPosition();
        radioactiveParticleOrientation = radData -> SetupOrientation();
        radioactiveParticleTime = radData -> GetInitialTime() * s;

        fParticleGun -> SetParticlePosition(radioactiveParticlePosition);
        fParticleGun -> SetParticleMomentumDirection(radioactiveParticleOrientation);
        fParticleGun -> SetParticleEnergy(radioactiveParticleEnergy);
        fParticleGun -> SetParticleDefinition(radioactiveParticle);
        fParticleGun -> SetParticleCharge(radioactiveParticleCharge);
        fParticleGun -> SetParticleTime(radioactiveParticleTime);
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
