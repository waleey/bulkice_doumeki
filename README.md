!!Welcome to bulkice_doumeki simulation!! *Read Me is not up to date!!*

This is a GEANT4 based simulation of different optical modules currently in use and to be deployed in future in IceCube Neutrino Observatory in South Pole. Currently available optical modules for simulation are MDOM, LOM16, LOM18, PDOM, D-Egg, WOM (currently under development). Except for WOM, the optical module simulation was initially written by M. Unland and C. Lozano, and later they were modified by me. The simulation also contains detailed depth dependent ice properties under antarctic ice sheet and it can simulation 40*40*40 cubic meter of ice at different temperature. Currently, it can simulate positron and electron flux from CCSN neutrinos, background radioactivity inside pressure vessel of MDOM, LOMs, and D-Egg, and photon wave with different Zenith Angle. It also accepts SNEWPY neutrino flux models in a python program called "merger" and use sntools to simulate the positron and electron flux from ibd and enees interactions. 
*Quick Start Guide*

-Make sure you have Geant4 installed in your local machine.  

-use "git clone https://github.com/waleey/bulkice_doumeki.git" to clone the rep on your local machine.  

-In build dir, change the variables inside the shell scripts according to your machine.  

-Modify CMakeLists.txt to make sure everything is added to the path correctly.   

-Go to bulkice_doumeki.cc and:  
	-change external variable "ghitsfilename" to where you wanna save the output file. Please follow the existing  
		format to avoid error.  
	-change gQEFile to the path of the QE file for the PMT on your local computer.  
	-change line 153 in OMSimDetectorConstruction.cc to the path of your build folder inside the bulkice_doumeki dir	
 -go to the build dir and type "cmake .."  
 
-then type "make -jN" inside build dir. N is the number of cores on your local machine. It helps compile faster.	
-follow instructions in bulkice_doumeki.cc to use different cmd line arguments to produce different  
outputs.


*NB: * I am still in the process of building and debugging it. So, please feel free to share any suggestion, criticism, or thoughts. My email: wkarim@u.rochester.edu slack: Waly M Z Karim  
