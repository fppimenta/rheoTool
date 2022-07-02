volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
    + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
);

MRF.makeRelative(phiHbyA);

if (p_rgh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustCorrPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

surfaceScalarField phig(thermo.hasGravity() ? -rAUf*thermo.ghf()*fvc::snGrad(thermo.rhok())*mesh.magSf() : phiHbyA*0.);

phiHbyA += phig;

tmp<volScalarField> rAtU(rAU);
 
if (simplec)
{
  rAtU = 1.0/(1.0/rAU - UEqn.H1());
  phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p_rgh)*mesh.magSf();
      
  if (thermo.hasGravity())
    HbyA -= (rAU - rAtU())*fvc::reconstruct(fvc::snGrad(p_rgh)*mesh.magSf());
  else
    HbyA -= (rAU - rAtU())*fvc::grad(p_rgh);
        
  // Update the pressure BCs to ensure flux consistency
  constrainPressure(p_rgh, U, phiHbyA, rAtU(), MRF);
}
else
{
  // Update the pressure BCs to ensure flux consistency
  constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);  
}
 
tUEqn.clear();
 
// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
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
   
    spSolverp_rghFluid[i].solve(p_rghEqn); 
                     
    if (pimple.finalNonOrthogonalIter())
    {
       phi = phiHbyA - p_rghEqn.flux();
     
       // Explicitly relax pressure for momentum corrector
       p_rgh.relax();

       if (thermo.hasGravity())
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
    }
}

#include "continuityErrs.H"
 
U.correctBoundaryConditions();
fvConstraints.constrain(U);
 
if (thermo.hasGravity())
{
  p == p_rgh + thermo.rhok()*thermo.gh();

  if (p_rgh.needReference())
  {
    p += dimensionedScalar
    (
      "p",
       p.dimensions(),
       pressureReference.refValue()
     - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - thermo.rhok()*thermo.gh();
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
