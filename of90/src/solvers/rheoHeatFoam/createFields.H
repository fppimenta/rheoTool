#include "createRDeltaT.H"

autoPtr<fluidThermoModel> thermo(fluidThermoModel::New(word("thermo"+mesh.name()), mesh));

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
 
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo->hasGravity()
    ?
    (p_rgh + thermo->rhok()*thermo->gh())()
    :
    p_rgh
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

pressureReference pressureReference(p, p_rgh, mesh.solutionDict().subDict("SIMPLE"));

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - thermo->rhok()*thermo->gh();
}

mesh.setFluxRequired(p_rgh.name());

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


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
   
