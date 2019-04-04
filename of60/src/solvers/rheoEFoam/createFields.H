Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

#include "createMRF.H"
#include "createFvOptions.H"

// Create matrix solver
#include "createSolver.H" 

// Create constitutive equation
constitutiveModel constEq(U, phi);

// Create the electric model    
EDFModel elecM(phi);
 
IOdictionary cttProperties
(
    IOobject
    (
        "constitutiveProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    )
);
   
bool   sPS = cttProperties.subDict("passiveScalarProperties").lookupOrDefault<Switch>("solvePassiveScalar", false);

autoPtr<volScalarField>  C 
(
  sPS != true
  ? 
  NULL
  :
  new volScalarField
  (
    IOobject
    (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
 )
);

   
