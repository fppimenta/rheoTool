#include "createRDeltaT.H"
 
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
 
// Create matrix solver
#include "createSolver.H" 

// Create constitutive equation
filmModel film(U, phi);
 
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

dimensionedVector g_(cttProperties.lookup("g"));
 
