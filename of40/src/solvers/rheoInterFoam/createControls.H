bool correctPhi
(
    pimple.dict().lookupOrDefault<Switch>("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
);

// Read extra-controls

int    nInIter = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<int>("nInIter", 1);

bool   simplec = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<Switch>("SIMPLEC", true);
 
