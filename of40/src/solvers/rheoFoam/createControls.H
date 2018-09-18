simpleControl simple(mesh);
    
#include "createTimeControls.H"

bool correctPhi
(
    simple.dict().lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    simple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

// Read extra-controls
int    nInIter = simple.dict().lookupOrDefault<int>("nInIter", 1);
  
