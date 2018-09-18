{
    volScalarField rAU = 1.0/UEqn.A();
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    U = rAU*UEqn.H();
    surfaceScalarField phiU
    (
       "phiU",
       mesh.changing() == true ?
       (fvc::interpolate(U) & mesh.Sf())
       :
       (fvc::interpolate(U) & mesh.Sf())
     + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    if (pd.needReference())
    {
        adjustPhi(phi, U, pd);
    }

    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );
    
    phi = phiU + phig;
        
    tmp<volScalarField> rAtU(rAU);

    if (simplec)
    {
      rAtU = 1.0/(1.0/rAU - UEqn.H1());
      rAtU().rename("rAtU");
    
      phi += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(pd)*mesh.magSf();
      U -= (rAU - rAtU())*fvc::grad(pd);
    }


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAtU(), pd, "laplacian(pd|(ap-H1))") == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);
        
        if (simplec)
        {
          pdEqn.solve();
        }
        else
        {
          pdEqn.solve(mesh.solutionDict().solver(pd.select(pimple.finalInnerIter())));
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }
    }

    if (simplec)
     {
          U += rAU*fvc::reconstruct(phig/rAUf) - rAtU()*fvc::grad(pd);
     }
    else
     {
          U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
     }
     
    U.correctBoundaryConditions();

    #include "continuityErrs.H"

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
