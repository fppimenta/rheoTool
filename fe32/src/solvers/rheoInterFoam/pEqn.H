{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    adjustPhi(phiU, U, pd);
        
    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf()
    );
    
    phi = phiU + phig;
        
    tmp<volScalarField> rAtU(rUA);

    rAtU = 1.0/(1.0/rUA - UEqn.H1());
    phi += fvc::interpolate(rAtU() - rUA)*fvc::snGrad(pd)*mesh.magSf();
    U -= (rUA - rAtU())*fvc::grad(pd);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAtU(), pd, "laplacian(pd|(ap-H1))") == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct(phig/rUAf) - rAtU()*fvc::grad(pd);
    U.correctBoundaryConditions();
}
