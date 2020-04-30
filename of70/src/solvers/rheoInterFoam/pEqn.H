{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
       fvc::makeRelative(phiHbyA, U);
       adjustCorrPhi(phiHbyA, U, p_rgh);
       fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
      (
          mixture.surfaceTensionForce()
        - ghf*fvc::snGrad(rho)
      )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;
    
    tmp<volScalarField> rAtU(rAU);

    if (simplec)
    {
      rAtU = 1.0/(1.0/rAU() - UEqn.H1());
      phiHbyA += fvc::interpolate(rAtU() - rAU())*fvc::snGrad(p_rgh)*mesh.magSf();
      HbyA -= (rAU() - rAtU()) * fvc::reconstruct( fvc::snGrad(p_rgh) * mesh.magSf() );
      
      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p_rgh, U, phiHbyA, rAtU(), MRF);
    }
    else
    {
      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);  
    }
 
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAtU(), p_rgh, "laplacian(p|rAtU)") == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        spSolverprgh->solve(p_rghEqn); 

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            if (simplec)
            {
               surfaceScalarField rAtUf("rAtUf", fvc::interpolate(rAtU()));
               U = HbyA + rAU()*fvc::reconstruct(phig/rAUf) - rAtU()*fvc::reconstruct(p_rghEqn.flux()/rAtUf);
            }
            else
            {
               U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            }
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
