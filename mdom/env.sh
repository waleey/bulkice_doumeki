#############################################
# 
# Geant4 environment parameters
# you must have correct path for following parameters.
# you may need to modify #data section too.
#

export G4BUILD=/Users/walu/icecube/my_geant4/geant4-v11.0.0-install
if [ "$G4BUILD" != "/Users/walu/icecube/my_geant4/geant4-v11.0.0-install" ]; then
    echo "Error: G4BUILD is not set to the correct path."
    exit 1  # Terminate the script with an error code
fi #terminate the run if the installation directory isn't set up properly!
#export G4BUILD=/cvmfs/icecube.opensciencegrid.org/users/hoshina/RHEL7/geant4.10.06.p02_mod_build
#export G4DATA=$G4BUILD/data/
#export G4VIS_BUILD_RAYTRACER_DRIVER=1
#export G4VIS_USE_RAYTRACER=1
source $G4BUILD/share/Geant4/geant4make/geant4make.sh
source $G4BUILD/bin/geant4.sh

#Only uncomment these lies if on mac:
OS=`uname -o`
if [[ "${OS}" == "Darwin" ]]; then
    export PKG_CONFIG_PATH=/opt/homebrew/opt/qt@5/lib/pkgconfig/
    export PATH=/opt/homebrew/opt/qt@5/bin:$PATH
fi
