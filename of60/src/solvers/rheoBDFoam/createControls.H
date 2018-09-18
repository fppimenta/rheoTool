simpleControl simple(mesh);
    
#include "createTimeControls.H"

bool correctPhi
(
    simple.dict().lookupOrDefault("correctPhi", mesh.dynamic())
);

bool checkMeshCourantNo
(
    simple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
 (
     simple.dict().lookupOrDefault("moveMeshOuterCorrectors", false)
 );

// Read extra-controls
int    nInIter = simple.dict().lookupOrDefault<int>("nInIter", 1);

bool   solveFluid = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<Switch>("solveFluid", true); 

bool   solveElecM = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<Switch>("solveElecM", false);
