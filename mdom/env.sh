#############################################
# 
# Geant4 environment parameters
# you must have correct path for following parameters.
# you may need to modify #data section too.
#

export G4BUILD=/home/waly/software/geant4-v11.1.2-install
#export G4BUILD=/cvmfs/icecube.opensciencegrid.org/users/hoshina/RHEL7/geant4.10.06.p02_mod_build
#export G4DATA=$G4BUILD/data/
#export G4VIS_BUILD_RAYTRACER_DRIVER=1
#export G4VIS_USE_RAYTRACER=1
source $G4BUILD/share/Geant4/geant4make/geant4make.sh
source $G4BUILD/bin/geant4.sh

