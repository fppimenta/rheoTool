#include "createTimeControls.H"

bool correctPhi
(
    simple.dict().lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    simple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

int    nInIter = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<int>("nInIter", 1);
 
