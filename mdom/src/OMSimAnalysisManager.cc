#include "OMSimAnalysisManager.hh"
#include "G4ios.hh"
//since Geant4.10: include units manually
#include "G4SystemOfUnits.hh"

extern G4int gDOM;
extern G4String ghitsfilename;

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
            datafile << stats_event_id.at(i) << "\t";
            datafile << std::fixed << stats_hit_time.at(i) /ns << "\t";
            //datafile << stats_photon_flight_time.at(i) /ns << "\t";
            //datafile << stats_photon_track_length.at(i) << "\t";
            datafile << stats_photon_energy.at(i) << "\t";
            datafile << stats_PMT_hit.at(i) << "\t";
            //datafile << stats_event_distance.at(i) << "\t";
            datafile << stats_photon_position.at(i).x()/m << "\t";
            datafile << stats_photon_position.at(i).y()/m << "\t";
            datafile << stats_photon_position.at(i).z()/m << "\t";
            datafile << stats_vertex_position.at(i).x()/m << "\t";
            datafile << stats_vertex_position.at(i).y()/m << "\t";
            datafile << stats_vertex_position.at(i).z()/m << "\t";
            datafile << stats_positron_id.at(i) << "\t";
            datafile << stats_survived_qe.at(i) << "\t";
           // datafile << stats_photon_direction.at(i).x() << "\t";
            //datafile << stats_photon_direction.at(i).y() << "\t";
            //datafile << stats_photon_direction.at(i).z() << "\t";
            //datafile << stats_photon_position.at(i).mag() / m ;
            datafile << stats_creator.at(i) << "\t";
            datafile << G4endl;
        }
	}
	else
	{
        G4cout << "********Failed to open " << ghitsfilename << " file*******" << G4endl;
	}


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

}

void OMSimAnalysisManager::SetRunId(G4int RunId)
{
    fRunId = RunId;
}
