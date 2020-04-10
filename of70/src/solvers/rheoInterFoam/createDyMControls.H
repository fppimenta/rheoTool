#include "createControl.H"
#include "createTimeControls.H"

bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", mesh.dynamic())
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault("moveMeshOuterCorrectors", false)
);

// Read extra-controls

int    nInIter = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<int>("nInIter", 1);

bool   simplec = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<Switch>("SIMPLEC", true);
 
