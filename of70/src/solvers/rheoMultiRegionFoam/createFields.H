#include "createFluidFields.H"
#include "createSolidFields.H"

//- Coupled T 
IOdictionary fvSolutionGen
(
  IOobject
  (
    "fvSolution",
    runTime.system(),
    runTime,
    IOobject::MUST_READ,
    IOobject::NO_WRITE,
    false
  )
); 

bool solveCoupledT = fvSolutionGen.subDict("coupledSolvers").subDict("T").lookupOrDefault<Switch>("solveCoupledT", true);

// Set flag at coupled T BC. The coeffs handling depends on that
if (solveCoupledT)
{
  //... for fluid
  forAll(fluidRegions, i)
  {
    volScalarField& T = const_cast<volScalarField&>(thermoFluid[i].T());
    forAll(T.boundaryField(), pI)
    {
      if (isType<coupledTFvPatchScalarField>(T.boundaryField()[pI])) 
        const_cast<coupledTFvPatchScalarField&>
        (
          refCast<const coupledTFvPatchScalarField>
          ( 
            T.boundaryField()[pI]
          )
        ).setFlagCoupled(true); 
    }  
    
    // Initialize kappa at interface
    T.correctBoundaryConditions();
  }
 
  //... for solid
  forAll(solidRegions, i)
  {  
    volScalarField& T = const_cast<volScalarField&>(thermoSolid[i].T());
    forAll(T.boundaryField(), pI)
    {
      if (isType<coupledTFvPatchScalarField>(T.boundaryField()[pI])) 
       const_cast<coupledTFvPatchScalarField&>
       (
          refCast<const coupledTFvPatchScalarField>
          ( 
            T.boundaryField()[pI]
          )
       ).setFlagCoupled(true);   
    }
    
    // Initialize kappa at interface
    T.correctBoundaryConditions();
  } 
}


//- Create coupled solver if needed
autoPtr<coupledSolver>  cpsT
(
  solveCoupledT
  ? 
  new coupledSolver("T","multiRegionMesh",runTime,fvSolutionGen)
  :
  nullptr
);

if (solveCoupledT)
{
  // Insert mesh
  //... fluid
  forAll(fluidRegions, i)
  {
    cpsT().insertMesh(fluidRegions[i]);
  }
 
  //... for solid
  forAll(solidRegions, i)
  {
    cpsT().insertMesh(solidRegions[i]);
  }
 
  // Insert T fields
  //... for fluid
  forAll(fluidRegions, i)
  {
    volScalarField& T = const_cast<volScalarField&>(thermoFluid[i].T());
    cpsT().insertField(T);
  }
 
  //... for solid
  forAll(solidRegions, i)
  {  
    volScalarField& T = const_cast<volScalarField&>(thermoSolid[i].T());
    cpsT().insertField(T);
  } 
}
    
 
 
