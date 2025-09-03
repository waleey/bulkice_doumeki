#include "OMSimAnalysisManager.hh"
#include "G4ios.hh"
//since Geant4.10: include units manually
#include "G4SystemOfUnits.hh"

extern G4int gDOM;
extern G4String ghitsfilename;
extern G4double gAngle;
extern G4double gPhotonSim;
extern G4int gVesselCount;
extern G4int gWLSCount;
extern G4int gPMTBodyCount;
extern G4int gWLSAndHitCount;

OMSimAnalysisManager::OMSimAnalysisManager(){
    std::cerr << "OMSimAnalysisManager is generated" << std::endl;}

OMSimAnalysisManager::~OMSimAnalysisManager(){}

void OMSimAnalysisManager::Write()
{
	if(datafile.is_open())
	{
        G4cout << "+++++++++++++The size of this event is " << stats_event_id.size() << " ++++++++++++" << G4endl;
        //datafile << "Run ID: " << fRunId << G4endl;
        //datafile << "Event Size: " << stats_event_id.size() << G4endl;
        for (int i = 0; i < (int) stats_event_id.size(); i++)
        {
            //datafile << stats_event_id.at(i) << "\t";
            datafile << fRunId << "\t";
            datafile << std::fixed << stats_hit_time.at(i) /ns << "\t";
            //datafile << stats_photon_flight_time.at(i) /ns << "\t";
            //datafile << stats_photon_track_length.at(i) << "\t";
            datafile << stats_photon_energy.at(i) << "\t";
            datafile << stats_PMT_hit.at(i) << "\t";
			datafile << stats_dom_id.at(i) << "\t"; //for multiple DOM simulation	
            //datafile << stats_event_distance.at(i) << "\t";
            datafile << stats_photon_position.at(i).x()/m << "\t";
            datafile << stats_photon_position.at(i).y()/m << "\t";
            datafile << stats_photon_position.at(i).z()/m << "\t";
            datafile << stats_vertex_position.at(i).x()/m << "\t";
            datafile << stats_vertex_position.at(i).y()/m << "\t";
            datafile << stats_vertex_position.at(i).z()/m << "\t";
            datafile << stats_positron_id.at(i) << "\t";
            datafile << stats_survived_qe.at(i) << "\t";
            if(gPhotonSim)
            {
                datafile << gAngle << "\t";
            }
           // datafile << stats_photon_direction.at(i).x() << "\t";
            //datafile << stats_photon_direction.at(i).y() << "\t";
            //datafile << stats_photon_direction.at(i).z() << "\t";
            //datafile << stats_photon_position.at(i).mag() / m ;
            //datafile << stats_creator.at(i) << "\t";
            datafile << G4endl;
        }
	}
	else
	{
        G4cout << "********Failed to open " << ghitsfilename << " file*******" << G4endl;
	}

//	std::cout << "Photon reaching Pressure Vessel: " << gVesselCount << std::endl
//	<< "Photon getting absorbed by PMT no_optic: " << gPMTBodyCount << std::endl
//	<< "Photon getting absorbed by WLS material: " << gWLSCount << std::endl
//	<< "Photon making hit and produced from WLS: " << gWLSAndHitCount << std::endl;

	gWLSAndHitCount = 0;
	gVesselCount = 0;
	gPMTBodyCount = 0;
	gWLSCount = 0;
}

void OMSimAnalysisManager::WriteAccept()
{
	int num_pmts;
	if (gDOM==0){num_pmts = 1;} //single PMT
	else if (gDOM==1){num_pmts = 24;} //mDOM
	else if (gDOM==2){num_pmts = 1;} //PDOM
	else if (gDOM==3){num_pmts = 16;} //LOM16
	else if (gDOM==4){num_pmts = 18;} //LOM18
	else if (gDOM==5){num_pmts = 2;} //DEGG
        else{num_pmts = 99;} //custom




	//int	pmthits[num_pmts+1] = {0};
        std::vector<int> pmthits(num_pmts+1, 0);
	int sum = 0;

	// repacking hits:
	for (int i = 0; i < (int) stats_PMT_hit.size(); i++) {
		pmthits[stats_PMT_hit.at(i)] += 1;
	}
	// wrinting collective hits
	for (int j = 0; j < num_pmts; j++) {
		datafile << "\t" << pmthits[j];
		sum += pmthits[j];
		pmthits[j] = 0;
	}
	datafile << "\t" << sum;
	datafile << G4endl;
}


void OMSimAnalysisManager::Reset()
{
	stats_event_id.clear();
	stats_photon_flight_time.clear();
	stats_photon_track_length.clear();
	stats_hit_time.clear();
	stats_photon_energy.clear();
	stats_PMT_hit.clear();
	stats_photon_direction.clear();
	stats_photon_position.clear();
	stats_vertex_position.clear();
	stats_positron_id.clear();
	stats_survived_qe.clear();
	stats_creator.clear();
	stats_event_distance.clear();
}

void OMSimAnalysisManager::SetRunId(G4int RunId)
{
    fRunId = RunId;
}
