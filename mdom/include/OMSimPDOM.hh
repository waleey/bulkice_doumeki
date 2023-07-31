#ifndef OMSimPDOM_h
#define OMSimPDOM_h 1
#include "abcDetectorComponent.hh"
#include "G4Orb.hh"

class pDOM : public abcDetectorComponent
    {
    public:
        pDOM(OMSimInputData* pData, G4bool pPlaceHarness = true);
        void Construction();
        G4bool mPlaceHarness;

        inline G4Orb* GetOuterSolid() { return lGlassSphereSolid; }
        inline G4Orb* GetInnerSolid() { return lGelSphereSolid; }


    private:
        OMSimPMTConstruction *mPMTManager;
        G4Orb* lGlassSphereSolid;
        G4Orb* lGelSphereSolid;

    };

#endif
