volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi, Uf))
);

MRF.makeRelative(phiHbyA);

if (p_rgh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustCorrPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

surfaceScalarField phig(thermo->hasGravity() ? -rAUf*thermo->ghf()*fvc::snGrad(thermo->rhok())*mesh.magSf() : phiHbyA*0.);

phiHbyA += phig;

tmp<volScalarField> rAtU(rAU);
 
if (simplec)
{
  rAtU = 1.0/(1.0/rAU - UEqn.H1());
  phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p_rgh)*mesh.magSf();
      
  if (thermo->hasGravity())
    HbyA -= (rAU - rAtU())*fvc::reconstruct(fvc::snGrad(p_rgh)*mesh.magSf());
  else
    HbyA -= (rAU - rAtU())*fvc::grad(p_rgh);
        
  // Update the pressure BCs to ensure flux consistency
  constrainPressure(p_rgh, U, phiHbyA, rAtU(), MRF);
  
  tUEqn.clear();
}
else
{
  // Update the pressure BCs to ensure flux consistency
  constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);  
}

// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
      fvm::laplacian(rAtU(), p_rgh, "laplacian(p|(ap-H1))") == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference
    (
            pressureReference.refCell(),
            getRefCellValue(p_rgh, pressureReference.refCell())
    );
   
    spSolverp_rgh->solve(p_rghEqn); 
                     
    if (simple.finalNonOrthogonalIter())
    {
       phi = phiHbyA - p_rghEqn.flux();
     
       // Explicitly relax pressure for momentum corrector
       p_rgh.relax();

       if (thermo->hasGravity())
       {
         if (simplec)
         {
           surfaceScalarField rAtUf("rAtUf", fvc::interpolate(rAtU()));
           U = HbyA + rAU*fvc::reconstruct(phig/rAUf) - rAtU()*fvc::reconstruct(p_rghEqn.flux()/rAtUf);
         }
         else
         {
           U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
         }
       }
       else
       {
         U = HbyA - rAtU()*fvc::grad(p_rgh);
       }
       
       U.correctBoundaryConditions();
       fvConstraints.constrain(U);
    }
}

#include "continuityErrs.H"

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

if (thermo->hasGravity())
{
  p == p_rgh + thermo->rhok()*thermo->gh();

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
}
else
{
  p == p_rgh;

  if (p_rgh.needReference())
  {
    p += dimensionedScalar
    (
      "p",
       p.dimensions(),
       pressureReference.refValue()
     - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p;
  }
}
