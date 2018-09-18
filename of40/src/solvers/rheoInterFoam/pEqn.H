{
    rAU = 1.0/UEqn.A();
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
      "phiHbyA",
      mesh.changing() == true ?
      fvc::flux(HbyA)
    + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, Uf())
      :
      fvc::flux(HbyA)
    + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
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
      rAtU = 1.0/(1.0/rAU - UEqn.H1());
      phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p_rgh)*mesh.magSf();
      HbyA -= (rAU - rAtU())*fvc::grad(p_rgh);
      
      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p_rgh, U, phiHbyA, rAtU());
    }
    else
    {
      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p_rgh, U, phiHbyA, rAUf);  
    }
  
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAtU(), p_rgh,"laplacian(p_rgh|(ap-H1))") == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
 
        if (simplec)
        {
          p_rghEqn.solve();
        }
        else
        {
          p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
        }
        
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            if (simplec)
             {
               U = HbyA + rAU*fvc::reconstruct(phig/rAUf) - rAtU()*fvc::grad(p_rgh);
             }
            else
             {
               U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
             }
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"
    
    if (mesh.changing())
    {
        Uf() = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf() += n*(phi/mesh.magSf() - (n & Uf()));
    }

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
}
