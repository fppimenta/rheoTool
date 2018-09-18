#include "createTimeControls.H"

bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

int    nInIter = pimple.dict().lookupOrDefault<int>("nInIter", 1);

bool   simplec = pimple.dict().lookupOrDefault<Switch>("SIMPLEC", true);
 
