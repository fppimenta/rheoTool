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

bool   simplec = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<Switch>("SIMPLEC", true);

int    nCorrectors = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<int>("nCorrectors", 1);

if (simplec)
 nCorrectors = 1;
 
