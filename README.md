!!Welcome to bulkice_doumeki simulation!!

The Optical Modules used in this simulation are originally written by M. Unland, C. Lozano, and L. Classen. The optical module design has been ported and modified to use with the bulkice simulation of 40*40*40 cubic meter of ice. So far, Supernova signal event can be simulated and studied with all the five available OM models, namely: MDOM, LOM18, LOM16, DOM, and single PMT. For the study of radioactive decay, only MDOM can be used but other modules will be added soon. 

*Quick Start Guide*

-Make sure you have Geant4 installed in your local machine. 
-use "git clone https://github.com/waleey/bulkice_doumeki.git" to clone the rep on your local machine.
-In build dir, change the variables inside the shell scripts according to your machine.
-Modify CMakeLists.txt to make sure everything is added to the path correctly. 
-Go to bulkice_doumeki.cc and:
	-change external variable "ghitsfilename" to where you wanna save the output file. Please follow the existing
		format to avoid error.
	-change gQEFile to the path of the QE file for the PMT on your local computer.
	-change line 153 in OMSimDetectorConstruction.cc to the path of your build folder inside the bulkice_doumeki dir	-go to the build dir and type "cmake .."
	-then type "make -jN" inside build dir. N is the number of cores on your local machine. It helps compile faster.	-follow instructions in bulkice_doumeki.cc to use different cmd line arguments to produce different
		outputs.


*NB: * I am still in the process of building and debugging it. So, please feel free to share any suggestion, criticism, or thoughts. My email: wkarim@u.rochester.edu slack: Waly M Z Karim  
